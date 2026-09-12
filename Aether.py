import torch
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
import numpy as np
import random
import os
import time
import matplotlib.pyplot as plt

# =================================================================
# 1. SEED
# =================================================================
def set_seed(seed):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# =================================================================
# 2. G-INIT (Геометрична инициализация)
# =================================================================
class GInit:
    def __init__(self, alpha=0.1):
        self.alpha = alpha

    def initialize_linear(self, layer):
        out_f, in_f = layer.weight.shape
        n = min(out_f, in_f)
        W = np.eye(out_f, in_f)
        for i in range(n - 1):
            W[i, i+1] = 0.1
            W[i+1, i] = 0.1
        W /= (np.linalg.norm(W) + 1e-8)
        return torch.FloatTensor(W * self.alpha)

def apply_ginit(model, alpha=0.1):
    g = GInit(alpha=alpha)
    for m in model.modules():
        if isinstance(m, nn.Linear):
            with torch.no_grad():
                m.weight.data.copy_(g.initialize_linear(m))
                if m.bias is not None:
                    m.bias.zero_()

# =================================================================
# 3. AETHER OPTIMIZER V4.1
# =================================================================
class AetherOptimizerV4(optim.Optimizer):
    T_AMBIENT    = 25.0
    T_TARGET     = 55.0
    T_MAX        = 85.0
    HEAT_COEFF   = 1.2
    DISSIPATION  = 0.04
    SENSITIVITY  = 0.3
    WARMUP_EP    = 5

    def __init__(self, params, lr=0.001, sectors=10):
        defaults = dict(lr=lr, sectors=sectors, alpha=0.08, c=0.8)
        super().__init__(params, defaults)
        self.global_step   = 0
        self.thermal_state = self.T_AMBIENT
        self.sigma_history = []
        self.current_epoch = 0

    @torch.no_grad()
    def step(self, current_loss=None):
        self.global_step += 1

        if self.current_epoch < self.WARMUP_EP:
            sigma = 0.0
        else:
            sigma = torch.sigmoid(
                torch.tensor((self.thermal_state - self.T_TARGET) * self.SENSITIVITY)
            ).item()
        self.sigma_history.append(sigma)

        current_writes = 0
        max_writes     = 0

        for group in self.param_groups:
            num_sectors = group['sectors']
            idx         = self.global_step % num_sectors
            lr          = group['lr']
            alpha       = group['alpha']
            c           = group['c']

            for p in group['params']:
                if p.grad is None: 
                    continue
                st = self.state[p]

                if 'm' not in st:
                    st['m']        = torch.zeros_like(p)
                    st['v']        = torch.zeros_like(p)
                    st['last_upd'] = torch.full_like(
                        p, self.global_step - num_sectors, dtype=torch.long
                    )

                st['m'].mul_(0.9).add_(p.grad, alpha=0.1)
                st['v'].mul_(0.999).addcmul_(p.grad, p.grad, value=0.001)
                b1  = 1 - 0.9   ** self.global_step
                b2  = 1 - 0.999 ** self.global_step
                upd = (st['m'] / b1) / ((st['v'] / b2).sqrt() + 1e-8)

                if sigma < 0.99:
                    p.data.add_(upd, alpha=-lr * (1.0 - sigma))

                n     = p.numel()
                chunk = max(1, n // num_sectors)
                start = idx * chunk
                end   = min(start + chunk, n)

                delta_t = (self.global_step - st['last_upd'].view(-1)[start:end]).float().clamp(min=1)
                M = (1.0 + alpha * (c**2) * 5.0 / (delta_t.sqrt() + 1e-8)).clamp(max=3.0)

                p.data.view(-1)[start:end].add_(
                    upd.view(-1)[start:end] * M,
                    alpha=-lr * sigma
                )
                st['last_upd'].view(-1)[start:end] = self.global_step

                current_writes += int(n * (1.0 - sigma) + (end - start) * sigma)
                max_writes     += n

        sparsity_ratio = current_writes / (max_writes + 1e-8)
        self.thermal_state += (
            sparsity_ratio * self.HEAT_COEFF
            - self.DISSIPATION * (self.thermal_state - self.T_AMBIENT)
        )
        self.thermal_state = min(self.thermal_state, self.T_MAX)

        return self.thermal_state, (1.0 - sparsity_ratio) * 100.0, sigma

# =================================================================
# 4. МОДЕЛ И ДАННИ
# =================================================================
class SmallCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(3, 32, 3, padding=1), nn.ReLU(), nn.MaxPool2d(2),
            nn.Conv2d(32, 64, 3, padding=1), nn.ReLU(), nn.MaxPool2d(2),
            nn.Flatten(),
            nn.Linear(64*8*8, 256), nn.ReLU(),
            nn.Linear(256, 10)
        )
    def forward(self, x): 
        return self.net(x)

def get_loaders(device):
    t_train = transforms.Compose([
        transforms.RandomHorizontalFlip(),
        transforms.RandomCrop(32, padding=4),
        transforms.ToTensor(),
        transforms.Normalize((0.4914,0.4822,0.4465), (0.2023,0.1994,0.2010))
    ])
    t_test = transforms.Compose([
        transforms.ToTensor(),
        transforms.Normalize((0.4914,0.4822,0.4465), (0.2023,0.1994,0.2010))
    ])
    train = torch.utils.data.DataLoader(
        torchvision.datasets.CIFAR10('./data', train=True, download=True, transform=t_train),
        batch_size=128, shuffle=True, num_workers=2, pin_memory=True
    )
    test = torch.utils.data.DataLoader(
        torchvision.datasets.CIFAR10('./data', train=False, transform=t_test),
        batch_size=256, num_workers=2, pin_memory=True
    )
    return train, test

def evaluate(model, loader, device):
    model.eval()
    correct = 0
    with torch.no_grad():
        for x, y in loader:
            correct += model(x.to(device)).argmax(1).eq(y.to(device)).sum().item()
    return correct / 10000 * 100

# =================================================================
# 5. СТАТИСТИЧЕСКИ БЕНЧМАРК
# =================================================================
def run_one(seed, epochs, device, train_loader, test_loader, use_aether=True):
    set_seed(seed)
    model = SmallCNN().to(device)

    if use_aether:
        apply_ginit(model)
        opt = AetherOptimizerV4(model.parameters(), lr=0.001, sectors=10)
    else:
        opt = optim.Adam(model.parameters(), lr=0.001)

    criterion = nn.CrossEntropyLoss()
    best_acc  = 0.0

    for ep in range(epochs):
        model.train()
        if use_aether:
            opt.current_epoch = ep

        for x, y in train_loader:
            x, y = x.to(device, non_blocking=True), y.to(device, non_blocking=True)
            opt.zero_grad(set_to_none=True)
            loss = criterion(model(x), y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()

        acc = evaluate(model, test_loader, device)
        if acc > best_acc: 
            best_acc = acc

    return best_acc

def run_statistical_benchmark():
    device  = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    EPOCHS  = 40
    SEEDS   = [42, 123, 777, 1337, 2024]
    print(f"\nDevice: {device} | Epochs: {EPOCHS} | Seeds: {SEEDS}")
    print("СТАТИСТИЧЕСКА ВАЛИДАЦИЯ – 5 независими runs")
    print("=" * 65)

    train_loader, test_loader = get_loaders(device)

    adam_results   = []
    aether_results = []

    for i, seed in enumerate(SEEDS):
        print(f"\n── Run {i+1}/5 | Seed: {seed} ──────────────────────────")
        t0 = time.time()

        acc_adam = run_one(seed, EPOCHS, device, train_loader, test_loader, use_aether=False)
        adam_results.append(acc_adam)
        print(f"  Adam:   {acc_adam:.2f}%")

        acc_ae = run_one(seed, EPOCHS, device, train_loader, test_loader, use_aether=True)
        aether_results.append(acc_ae)
        print(f"  Aether: {acc_ae:.2f}%  ({acc_ae-acc_adam:+.2f}%)")
        print(f"  Време:  {(time.time()-t0)/60:.1f} мин")

    adam_arr   = np.array(adam_results)
    aether_arr = np.array(aether_results)
    diff_arr   = aether_arr - adam_arr

    adam_mean, adam_std     = adam_arr.mean(),   adam_arr.std()
    aether_mean, aether_std = aether_arr.mean(), aether_arr.std()
    diff_mean, diff_std     = diff_arr.mean(),   diff_arr.std()

    t_stat = diff_mean / (diff_std / np.sqrt(len(SEEDS)) + 1e-8)
    significant = abs(t_stat) > 2.776

    print("\n" + "=" * 65)
    print("  ФИНАЛНИ СТАТИСТИЧЕСКИ РЕЗУЛТАТИ")
    print("=" * 65)
    print(f"  Adam   mean±std: {adam_mean:.2f}% ± {adam_std:.2f}%")
    print(f"  Aether mean±std: {aether_mean:.2f}% ± {aether_std:.2f}%")
    print(f"  Разлика:         {diff_mean:+.2f}% ± {diff_std:.2f}%")
    print(f"  t-статистика:    {t_stat:.3f}")
    print(f"  Статистически значимо (p<0.05): {'✅ ДА' if significant else '❌ НЕ'}")
    print("=" * 65)

if __name__ == "__main__":
    run_statistical_benchmark()