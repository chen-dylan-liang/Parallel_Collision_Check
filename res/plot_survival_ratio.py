import pandas as pd
import matplotlib.pyplot as plt

file_path = 'benchmark_4_broad_effect.txt'
with open(file_path, 'r') as f:
    lines = f.readlines()

data = {'n': [], 'pairs': [], 'pairs_to_narrow': []}

for line in lines:
    line = line.strip()
    if line.startswith('n='):
        current_n = int(line.split('=')[1])
        data['n'].append(current_n)
    elif line.startswith('#pairs='):
        pairs = int(line.split('=')[1])
        data['pairs'].append(pairs)
    elif line.startswith('# pairs_to_narrow='):
        # Only record the first occurrence per n
        if len(data['pairs_to_narrow']) < len(data['n']):
            data['pairs_to_narrow'].append(int(line.split('=')[1]))

# Compute survival ratio
survival_ratio = [pn / p for pn, p in zip(data['pairs_to_narrow'], data['pairs'])]

# Create DataFrame
df = pd.DataFrame({'n': data['n'], 'survival_ratio': survival_ratio})

# Plot
plt.figure(figsize=(10, 6))
plt.plot(df['n'], df['survival_ratio'], color='blue', linewidth=3)
plt.xlabel('# Primitives', fontsize=14)
plt.ylabel('Survival Ratio', fontsize=14)
plt.title('Effect of Broad Phase', fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.savefig('survival_ratio.pdf', format='pdf', bbox_inches='tight')
