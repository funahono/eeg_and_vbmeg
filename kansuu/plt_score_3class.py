import numpy as np
import matplotlib.pyplot as plt

def visualize_decision_scores(decision: np.ndarray, test_labels: np.ndarray, label_list) -> plt.Figure:

    max_val = np.amax(decision)
    min_val = np.amin(decision)
    range_val = max(abs(max_val), abs(min_val)) + 0.5
    toumeido = 0.2
    colors = ['blue', 'green', 'red']

    decision_0 = decision[test_labels == 0]
    decision_1 = decision[test_labels == 1]
    decision_2 = decision[test_labels == 2]

    class_0 = label_list[0]
    class_1 = label_list[1]
    class_2 = label_list[2]

    fig, axes = plt.subplots(3, 2, figsize=(8, 6))  # 3行2列のサブプロット

    ## {class_0}
    y = np.arange(0, len(decision_0), 1)
    axes[0, 0].plot(decision_0[:, 0], y, 'o', markersize=4, color=colors[0])
    axes[0, 0].set_title(f'{class_1} vs {class_0}')
    axes[0, 0].set_xlim(-range_val, range_val)
    axes[0, 0].axvline(0, color='black', linestyle='--')
    axes[0, 0].set_yticks([])
    axes[0, 0].set_yticklabels([])
    axes[0, 0].set_xticks([])
    axes[0, 0].set_xticklabels([])
    axes[0, 0].axvspan(-range_val, 0, color=colors[1], alpha=toumeido)
    axes[0, 0].axvspan(0, range_val, color=colors[0], alpha=toumeido)
    axes[0, 0].spines['right'].set_visible(False)

    axes[0, 1].plot(-decision_0[:, 1], y, 'o', markersize=4, color=colors[0])
    axes[0, 1].set_title(f'{class_0} vs {class_2}')
    axes[0, 1].set_xlim(-range_val, range_val)
    axes[0, 1].axvline(0, color='black', linestyle='--')
    axes[0, 1].set_yticks([])
    axes[0, 1].set_yticklabels([])
    axes[0, 1].set_xticks([])
    axes[0, 1].set_xticklabels([])
    axes[0, 1].axvspan(-range_val, 0, color=colors[0], alpha=toumeido)
    axes[0, 1].axvspan(0, range_val, color=colors[2], alpha=toumeido)
    axes[0, 1].spines['left'].set_visible(False)

    ## {class_1}
    y = np.arange(0, len(decision_1), 1)
    axes[1, 0].plot(-decision_1[:, 0], y, 'o', markersize=4, color=colors[1])
    axes[1, 0].set_title(f'{class_0} vs {class_1}')
    axes[1, 0].set_xlim(-range_val, range_val)
    axes[1, 0].axvline(0, color='black', linestyle='--')
    axes[1, 0].set_yticks([])
    axes[1, 0].set_yticklabels([])
    axes[1, 0].set_xticks([])
    axes[1, 0].set_xticklabels([])
    axes[1, 0].axvspan(-range_val, 0, color=colors[0], alpha=toumeido)
    axes[1, 0].axvspan(0, range_val, color=colors[1], alpha=toumeido)
    axes[1, 0].spines['right'].set_visible(False)
    axes[1, 1].spines['left'].set_visible(False)

    axes[1, 1].plot(-decision_1[:, 2], y, 'o', markersize=4, color=colors[1])
    axes[1, 1].set_title(f'{class_1} vs {class_2}')
    axes[1, 1].set_xlim(-range_val, range_val)
    axes[1, 1].axvline(0, color='black', linestyle='--')
    axes[1, 1].set_yticks([])
    axes[1, 1].set_yticklabels([])
    axes[1, 1].set_xticks([])
    axes[1, 1].set_xticklabels([])
    axes[1, 1].axvspan(-range_val, 0, color=colors[1], alpha=toumeido)
    axes[1, 1].axvspan(0, range_val, color=colors[2], alpha=toumeido)

    ## {class_2}
    y = np.arange(0, len(decision_2), 1)
    axes[2, 0].plot(-decision_2[:, 1], y, 'o', markersize=4, color=colors[2])
    axes[2, 0].set_title(f'{class_0} vs {class_2}')
    axes[2, 0].set_xlim(-range_val, range_val)
    axes[2, 0].axvline(0, color='black', linestyle='--')
    axes[2, 0].set_yticks([])
    axes[2, 0].set_yticklabels([])
    axes[2, 0].set_xticks([])
    axes[2, 0].set_xticklabels([])
    axes[2, 0].axvspan(-range_val, 0, color=colors[0], alpha=toumeido)
    axes[2, 0].axvspan(0, range_val, color=colors[2], alpha=toumeido)
    axes[2, 0].spines['right'].set_visible(False)

    axes[2, 1].plot(decision_2[:, 2], y, 'o', markersize=4, color=colors[2])
    axes[2, 1].set_title(f'{class_2} vs {class_1}')
    axes[2, 1].set_xlim(-range_val, range_val)
    axes[2, 1].axvline(0, color='black', linestyle='--')
    axes[2, 1].set_yticks([])
    axes[2, 1].set_yticklabels([])
    axes[2, 1].set_xticks([])
    axes[2, 1].set_xticklabels([])
    axes[2, 1].axvspan(-range_val, 0, color=colors[2], alpha=toumeido)
    axes[2, 1].axvspan(0, range_val, color=colors[1], alpha=toumeido)
    axes[2, 1].spines['left'].set_visible(False)

    plt.subplots_adjust(wspace=0, hspace=0.3)

    return fig