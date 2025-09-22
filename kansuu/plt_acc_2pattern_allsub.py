import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import string

## plot ###############################################

def pltbar_average_allsub(data0, data1, label_list, chance_line, colors, calculate_average, legend):

    ## setting ############################################
    rcParams['font.family'] = 'DejaVu Sans'
    rcParams['font.size'] = 25

    width = 0.3

    chance_level = chance_line

    ## subject list #######################################
    alpha_list = list(string.ascii_uppercase) # list of subject name

    subject_list = []
    for i in range(len(data1)):
        subject_list.append(alpha_list[i])

    ## calculate average ##################################

    avg0 = sum(data0) / len(data0)
    avg1 = sum(data1) / len(data1)

    print('---------------------- average ----------------------')
    print(f'{label_list[0]}:{avg0}%')
    print(f'{label_list[1]}:{avg1}%')

    if calculate_average == 'true':
        data0 = np.append(data0, [0, avg0])
        data1 = np.append(data1, [0, avg1])
        subject_list = np.append(subject_list, ['','avg.'])

    ## plot setting ############################################

    figure = plt.figure(figsize=(12, 5))
    x = np.arange(len(subject_list))

    plt.ylabel('Accuracy (%)')
    plt.xlabel('Subject')

    ## plot result ##############################################

    x = np.arange(len(subject_list))
    width = 0.3 #line width

    for i in range(len(subject_list)):
        x0 = x[i] - width
        y0 = data0[i] - 0.5

        x1 = x[i]
        y1 = data1[i] - 0.5

        # plt.plot([x0, x1], [y0, y1], color='black', linestyle='--', linewidth=0.5)

    if legend == 'false':
        plt.bar(x - width/2, data0, width, color= colors[0], edgecolor = 'black')
        plt.bar(x + width/2, data1, width, color = colors[1], edgecolor = 'black')

    elif legend == 'true':
        plt.bar(x - width/2, data0, width, color= colors[0], edgecolor = 'black', label = label_list[0])
        plt.bar(x + width/2, data1, width, color = colors[1], edgecolor = 'black', label = label_list[1])
        plt.legend(loc = 'lower right')


    plt.xticks(x, subject_list)

    plt.ylim(0, 100)
    plt.yticks(np.arange(0, 101, 10))
    plt.subplots_adjust(right = 0.8)

    #ch level
    plt.axhline(y = chance_level, color = 'black', linestyle = '--', linewidth = 1.5)
    plt.close()

    return figure, avg0, avg1, subject_list

