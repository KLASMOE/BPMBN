import numpy as np
import matplotlib.pyplot as plt
import itertools


def parse_file(filename):
    data = []
    with open(filename, 'r') as fid:
        for line in fid:
            data_block = [list(map(float, fid.readline().strip().split()))
                          for _ in range(100)]

            info_array = np.array(list(map(float, line.split())))
            data_array = np.array(data_block)

            entry = {
                'info': info_array,
                'avg_triangulation_time': np.mean(data_array[:, 0]),
                'avg_computation_time': np.mean(data_array[:, 1]),
                'median_computation_time': np.median(data_array[:, 1]),
                'std_computation_time': np.std(data_array[:, 1])
            }
            data.append(entry)
    return data


def parse2(data, group_sizes):
    """
    按照指定的组大小分组，并计算每组的平均值。

    :param data: 由 parse_file 函数返回的数据列表。
    :param group_sizes: 一个列表，包含每个组的大小。
    :return: 一个列表，包含每个组的各统计值的平均值。
    """
    start_index = 0
    mean_results = []

    for size in group_sizes:
        chunk = data[start_index:start_index + size]
        start_index += size

        # 初始化用于累加各统计值的变量
        sum_avg_triangulation_time = 0
        sum_avg_computation_time = 0

        # 遍历当前块，并累加统计值
        for entry in chunk:
            sum_avg_triangulation_time += entry['avg_triangulation_time']
            sum_avg_computation_time += entry['avg_computation_time']

        # 计算当前块的每个统计值的平均，并添加到结果列表中
        n_items = len(chunk)
        mean_results.append({
            'mean_avg_triangulation_time':
                round(sum_avg_triangulation_time / n_items, 4),
            'mean_avg_computation_time':
                round(sum_avg_computation_time / n_items, 4),
        })

    return mean_results


def visualize(parsed_data, vertices, metric, metric_label, file_prefix,
              target_dir='./Images/'):
    markers = ['.', 'P', 'X', '*', 'o']
    linestyles = [(0, (3, 1, 1, 1, 1, 1)), ':', '-.', '--', '-']

    for vertex in vertices:
        markers_cycle = itertools.cycle(markers)
        linestyles_cycle = itertools.cycle(linestyles)

        plt.figure()
        lines_data = {}
        all_evidence_vars = set()

        for entry in parsed_data:
            if entry['info'][0] == vertex:
                ratio, evidence_vars = entry['info'][1], entry['info'][2]
                all_evidence_vars.add(evidence_vars)

                if ratio not in lines_data:
                    lines_data[ratio] = {'x': [], 'y': []}

                lines_data[ratio]['x'].append(evidence_vars)
                lines_data[ratio]['y'].append(entry[metric])

        for ratio, values in lines_data.items():
            plt.plot(values['x'], values['y'],
                     linestyle=next(linestyles_cycle),
                     marker=next(markers_cycle),
                     label=f"Cont. Vars. Prop. = {ratio:.2f}")

        plt.xticks(sorted(all_evidence_vars))
        plt.xlabel('Size of Evidence Set')
        plt.ylabel(metric_label)
        plt.legend(loc='best')
        plt.grid(False)

        eps_filename = f"{file_prefix}_{vertex}.eps"
        png_filename = f"{file_prefix}_{vertex}.png"
        plt.savefig(f'{target_dir}/{eps_filename}', format='eps', dpi=600)
        # plt.savefig(f'{target_dir}/{png_filename}', format='png', dpi=600)

        plt.show()
        plt.close()


if __name__ == "__main__":
    data_dir = './BPMBN/BPMBN/results/'
    data_file = data_dir + 'timecost.txt'
    parsed_data = parse_file(data_file)

    visualize(parsed_data, [50, 75, 100],
              'avg_computation_time', 'Average Computation Time (seconds)', 'avg')
    visualize(parsed_data, [50, 75, 100],
              'median_computation_time', 'Median Computation Time (seconds)', 'med')

    # Total average time costs
    # to find minimal strong triangulations and strong junction trees
    # for networks with 50, 75 and 100 variables
    group_sizes = [50, 75, 100]
    node_labels = ["50 nodes", "75 nodes", "100 nodes"]
    mean_results = parse2(parsed_data, group_sizes)

    for label, chunk_means in zip(node_labels, mean_results):
        print(f"{label}: {chunk_means}")

    # Time costs of network
    # with 100 variables 0.5 continuous variables proportion and 10 evidences
    no = 167
    print(parsed_data[no].get('info'))
    print('avg_computation_time: ',
          round(parsed_data[no].get('avg_computation_time'), 4))
    print('median_computation_time:',
          round(parsed_data[no].get('median_computation_time'), 4))
    print('std_computation_time: ',
          round(parsed_data[no].get('std_computation_time'), 4))
