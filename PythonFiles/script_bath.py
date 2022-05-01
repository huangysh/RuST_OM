# -*- coding: utf-8 -*-
# ======================================================================================================================
# Copyright 2021 School of Environmental Science and Engineering, Shanghai Jiao Tong University, Shanghai, PR. China.

# This script is used to calculated the cost of the rural sewage treatment system (contain the collection system and
# wastewater treatment system) of a considered area (a town, about 6000 households), when using facilities with various
# scales. And to obtain the variation curve between the investment (treatment system and collection system) and the
# number of wastewater treatment facilities.

# This script used a greedy algorithm to generate the optimal sewage collection system, and between two nodes we used
# the A-star algorithm to calculate the optimal pipeline. One of the most important target of this script is to obtain
# the most cost-effective solution for rural wastewater treatment. Therefore, the total pipeline cost (construction,
# operation and maintenance) was set as the weight between in A* algorithm, rather than distance.

# 用于批量处理数据, 如不同目标区域的模拟；同一目标区域不同覆盖率的模拟

# ----------------------------------------------------------------------------------------------------------------------
# GNU General Public License
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# ----------------------------------------------------------------------------------------------------------------------
# Author: ysh_Huang
# Contact: YuanshengHuang@hotmail.com
# Date: 22-05-2021
# Version: 0.1.1
# ======================================================================================================================

# general import
from __future__ import print_function

from functions import *
import time

import matplotlib.pyplot as plt

father_folder = r"E:\RuST\SSM\input_data_folder"
dem_folder = father_folder + '//' + 'DEM'
name_list_x = get_file(dem_folder, '.shp')

name_lis0 = ['Qian_B.shp', 'Zhe_A.shp', 'Xiang_D.shp', 'Liao_B.shp', 'Chuan_A.shp', 'Qian_D.shp', 'Min_A.shp',
             'E_D.shp', 'Zhe_D.shp', 'Qiong_B.shp', 'Hei_A.shp', 'JinSX_A.shp', 'JinSX_B.shp', 'Jin_A.shp', 'Su_B.shp',
             'JiHB_B.shp', 'YuHN_B.shp', 'Xiang_B.shp', 'Meng_A.shp', 'Xiang_C.shp', 'Xiang_C.shp', 'JiJL_B.shp']

name_list = ['Zhe_A.shp', 'Yue_D.shp', 'Qiong_B.shp', 'YuHN_B.shp']

# for c in [1, 0.98, 0.96, 0.94, 0.9, 0.85]:  # coverage
coverage_list = [1]   # user input

for c in coverage_list:  # coverage
    coverage = c  # 污水处理覆盖率

    # for name in name_list:
    for name in name_list:
        print('****************************************')
        print(name, end='')
        print(': ' + str(c), end='')
        print("  Time: ", end='')
        print(time.asctime(time.localtime(time.time())))
        print('****************************************')

        # set SSM parameters
        print("Setting SSM parameters...")

        script_folder = r'E:\RuST\SSM\RuST_simulation V0.1.1'  # input parameter 01
        input_data_folder = r"E:\RuST\SSM\input_data_folder"
        database_path = script_folder + "//" + "dataset"  # todo
        output_folder = r'E:\RuST\SSM\output_data_folder'

        # setting input data
        households = r'E:\RuST\RuST_points_datasets\clustered_households\shapeFile' + '//' + name
        DEM_points = r'E:\RuST\RuST_points_datasets\DEM100_P' + '//' + name
        ENVI_point = r'E:\RuST\RuST_points_datasets\envi_30' + '//' + name

        # setting output
        output_file = output_folder + "//" + str(c) + name + "result.xls"

        # setting builtin parameters folder
        standard_f = database_path + "//" + "standard.csv"
        technology_f = database_path + "//" + "technology.csv"
        price_onsite_f = database_path + "//" + "price_onsite.csv"
        pipe_type_f = database_path + "//" + "pipe_type.csv"
        pipe_para_f = database_path + "//" + "pipe_para.csv"
        price_f = database_path + "//" + "price.csv"
        cellsize = 100

        print("Reading SSM parameters files...")
        # readout database to a dictionary
        standard = read_file(standard_f, n=1)  # {id: [down, up, envi, level], ...}
        del standard['0']  # 删除表头
        price_onsite = read_file(price_onsite_f, n=1)  # {id: [down, up, cost], ...}  不含化粪池、隔油井
        pipe_type = read_file(pipe_type_f, n=1)  # {id: [diameter, price], ...}
        pipe_para = read_file(pipe_para_f, n=1)  # {id: [price], ...}
        price = read_file(price_f, n=1)  # {id: [price], [], ...}
        technology = read_file_list(technology_f, n=2)  # [[id, name, level, down, up, cost, O&M], ...]
        price_list = [price['earthwork'], pipe_type, price["drop"], price['inspection'], price['pump']]

        onsite_cost = []
        for p in price_onsite.values():
            a = [p[0], p[1], p[2] + price["septic"][0] + price["grease"][0]]
            onsite_cost.append(a)

        # get pipeline design parameters  # pipeline design parameters for sensitive analysis
        pop, charge = pipe_para["population"][0], pipe_para['sewage'][0]
        slope = pipe_para["slope"][0]
        depth_min, depth_max, wide = pipe_para["depth_min"][0], pipe_para["depth_max"][0], pipe_para["trench_wide"][0]
        sewer_ls, facility_ls = pipe_para["sewer_ls"][0], pipe_para["facility_ls"][0]

        print("Reading input data sets...")
        # reading input data
        households_list_o, _ = read_households_point(households, pop)
        households_list = coverage_rate(households_list_o, coverage)
        households_dict = list_to_dict(households_list)
        graph, _ = read_DEM_point(DEM_points, pop)
        ENVI, _ = read_ENVI_point(ENVI_point)

        print("--------------------------------------------------------------")
        print("Number of household: ", end='')
        print(len(households_list_o))
        print("Covered households:", end='')
        print(len(households_list))
        print("--------------------------------------------------------------")

        print("Ready to run Scenario Simulation Model... ")
        print(time.asctime(time.localtime(time.time())))

        final_out = []
        output_list_a, cluster_list = onsite_to_cluster(households_list, onsite_cost, technology, ENVI, standard,
                                                        facility_ls, sewer_ls, cellsize, pop, charge, n=0)

        for i in output_list_a:
            final_out.append(i)

        # data clean
        for i in cluster_list:
            id_list = list(set(i[7]))
            i[7] = id_list
        # assign elevation
        cluster_list = assign_height(cluster_list, graph, cellsize)

        print("--------------------------------------------------------------")
        print("Number of cluster: ", end='')
        print(len(cluster_list))
        print("--------------------------------------------------------------")

        print("Household merging preformed successfully!!!")
        print(time.asctime(time.localtime(time.time())))

        message = "There are %d cluster within the considered area." % len(cluster_list)
        print(message)

        print("Ready to merge cluster... ")
        # get sewer cost within cluster
        sewer_cost = [output_list_a[-1][5], output_list_a[-1][6]]

        print(time.asctime(time.localtime(time.time())))
        output_list_b = to_centralized(households_dict, cluster_list, graph, pipe_type, price_list, sewer_cost,
                                       technology,
                                       ENVI,
                                       standard, 1.5, depth_min, wide, depth_max, slope, sewer_ls, facility_ls, 50,
                                       cellsize, pop,
                                       charge, n=0)

        for j in output_list_b:
            final_out.append(j)

        print("Cluster merging preformed successfully!!!")
        print(time.asctime(time.localtime(time.time())))

        print('final_out:', end='')
        print(final_out)

        print("Ready to write out... ")
        output_dict = get_output_dict(final_out, len(households_list))

        write_to_file(output_dict, output_file)

        print("Write out successfully!!!")
        print("Ready to draw out... ")
        print(time.asctime(time.localtime(time.time())))

        # draw plot

        plt.title('SSM Results')
        plt.plot(output_dict["04_DD_new"], output_dict["12_Total_cost"], color="g", marker='o', label='Total Cost')
        plt.plot(output_dict["04_DD_new"], output_dict["10_Sewer_cost"], color="r", marker='^', label='Sewer Cost')
        plt.plot(output_dict["04_DD_new"], output_dict["11_WWTP_cost"], color="b", marker='*', label='Facility Cost')

        plt.legend()

        plt.xlabel('Degree of decentralization')
        plt.ylabel('Average Investment (CNY/Household)')
        plt.show()

        print(time.asctime(time.localtime(time.time())))
        print("SSM performed successfully")

        gc.collect()


if __name__ == '__main__':
    pass