import pandas as pd
import matplotlib.pyplot as plt

# Load and parse the data
file_path = 'benchmark_4.txt'

with open(file_path, 'r') as f:
    lines = f.readlines()


data = {
    'n': [],
    'serial_parry_narrow': [],
    'serial_our_narrow': [],
    'parallel_parry_narrow': [],
    'parallel_our_narrow': [],
    'serial_double': [],
    'parallel_double': [],
}


def time_to_ms(time_str):
    if 'µs' in time_str:
        return float(time_str.replace('µs', '')) / 1000
    elif 'ms' in time_str:
        return float(time_str.replace('ms', ''))
    elif 's' in time_str:
        return float(time_str.replace('s', '')) * 1000
    else:
        return float(time_str)


current_n = None
for line in lines:
    line = line.strip()
    if line.startswith('n='):
        current_n = int(line.split('=')[1])
        data['n'].append(current_n)
    else:
        for key in data.keys():
            if key != 'n' and line.startswith(key):
                time_str = line.split('=')[1]
                data[key].append(time_to_ms(time_str))


df = pd.DataFrame(data)


plt.figure(figsize=(10, 6))


style_map = {
    'parallel_double': {'color': 'blue', 'linestyle': '-', 'linewidth': 3},
    'serial_double': {'color': 'red', 'linestyle': '-', 'linewidth': 3},
    'parallel_our_narrow': {'color': 'green', 'linestyle': '-', 'linewidth': 3},
    'parallel_parry_narrow': {'color': 'purple', 'linestyle': '--', 'linewidth': 2},
    'serial_our_narrow': {'color': 'orange', 'linestyle': '--', 'linewidth': 2},
    'serial_parry_narrow': {'color': 'cyan', 'linestyle': '--', 'linewidth': 2},
}


for method in style_map.keys():
    plt.plot(df['n'], df[method],
             label=method,
             color=style_map[method]['color'],
             linestyle=style_map[method]['linestyle'],
             linewidth=style_map[method]['linewidth'])


plt.xlabel('# Primitives', fontsize=14)
plt.ylabel('Time (ms)', fontsize=14)
plt.title('Benchmark Comparisons (Logarithm Y-axis)', fontsize=16)
plt.grid(True)
yscale_style='log'
#plt.xscale('log')
plt.yscale(yscale_style)


sorted_methods = sorted(style_map.keys(), key=lambda m: df[m].iloc[-1], reverse=True)
idx_serial_double = sorted_methods.index('serial_double')
idx_parallel_our_narrow = sorted_methods.index('parallel_our_narrow')


handles, labels = plt.gca().get_legend_handles_labels()
sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: sorted_methods.index(x[1]))
handles, labels = zip(*sorted_handles_labels)

plt.legend(handles, labels, fontsize=12, loc='best')
plt.tight_layout()
plt.savefig("scale_"+yscale_style+".pdf", format='pdf', bbox_inches='tight')