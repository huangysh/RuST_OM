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

# This file contains all functions, which would be used in the ArcGIS.

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
import gc
import os
import math
import arcpy
import numpy as np
import pandas as pd

from copy import deepcopy


# ----------------------------------------------------------------------------------------------------------------------
# data readout
# ----------------------------------------------------------------------------------------------------------------------


def get_file(folder, marker=None):
    """
    This function is used to get the target file names within a certain folder.

    :param folder: the target folder path
    :param marker: special marker (file suffix) of file name with the target folder.

    :return:
    file_name_list: target file name within the target folder.
    """
    file_name_list = []
    readout = os.walk(folder)
    file_name = list(readout)[0][2]
    num = len(marker)
    for i in file_name:
        if marker == i[- num:]:
            file_name_list.append(i)
        else:
            continue
    return file_name_list


def read_households_point(shapefile, pop=3):
    """
    This function is used to readout point feature, and the attribute table of this point feature should contain fields
    such as X, Y, Z, Q, marker.

    :param shapefile: point feature, attribute table with fields X, Y, Z, Q, marker. 原始数据中的Q任然为人口，非流量。
    :param pop: 户均人口，默认为3

    :return:
    point_list: nested list, [[ID, X, Y, H, Q, M, 0.5, [ID], [envi, level], tech, c, o&m, []], [], ...] 0.5为最小埋深，米
    """
    point_list, rows, fields = [], arcpy.SearchCursor(shapefile), arcpy.ListFields(shapefile)
    spatial_ref = arcpy.Describe(shapefile).spatialReference
    ID, X, Y, Q, H, M = 100000, 0, 0, 0, 0, None
    id_x = None
    for row in rows:
        for field in fields:
            if field.name == "POINT_X" or field.name == "point_x":
                X = row.getValue(field.name)
            if field.name == "POINT_Y" or field.name == "point_y":
                Y = row.getValue(field.name)
            # if field.name == "Q" or field.name == "q":
            # Q = row.getValue(field.name)
            if field.name == "RASTERVALU" or field.name == "rastervalu" or field.name == "GRID_CODE" or field.name == "grid_code":
                H = row.getValue(field.name)
            if field.name == "markerL":
                M = row.getValue(field.name)
                id_x = "H%d" % ID  # 编号，起始点，H100000
        point_list.append([id_x, X, Y, H, pop, M, 0.5, [id_x], ['e', 'l'], None, 0, 0, [0, 0, 0, 0, 0]])
        ID += 1

    if len(point_list) < 1:
        raise Exception("EMPTY LIST: YOU GOT AN EMPTY LIST, PLEASE CHECK YOUR INPUT FILE!!!")
    return point_list, spatial_ref


def read_DEM_point(shapefile, pop):
    """
    Read spatial point feature and store as a list

    :param shapefile: point feature (shapefile) with fields， such as point_x, point_y and rastervalu
    :param pop: 户均人数

    :return:
    point_list: nested list, [[ID, X, Y, H, "Q", "M", 0.5, [ID], [envi, level], tech, c, o&m, []], [], ...] 0.5为最小埋深
    """
    point_list, rows, fields, ID = [], arcpy.SearchCursor(shapefile), arcpy.ListFields(shapefile), 10000000
    spatial_ref = arcpy.Describe(shapefile).spatialReference
    X, Y, H = 0, 0, 0
    ID_P = "D%d" % ID
    for row in rows:
        for field in fields:
            if field.name == "POINT_X" or field.name == "point_x":
                X = row.getValue(field.name)
            if field.name == "POINT_Y" or field.name == "point_y":
                Y = row.getValue(field.name)
            if field.name == "RASTERVALU" or field.name == "rastervalu" or field.name == "GRID_CODE" or field.name == "grid_code":
                H = row.getValue(field.name)
                ID_P = "D%d" % ID  # 编号，起始点，E10000000
        point_list.append([ID_P, X, Y, H, pop, "M", 0.5, [ID_P], ['e', 'l'], None, 0, 0, [0, 0, 0, 0, 0]])
        ID += 1

    if len(point_list) < 1:
        raise Exception("EMPTY LIST: YOU GOT AN EMPTY LIST, PLEASE CHECK YOUR INPUT FILE!!!")
    return point_list, spatial_ref


def read_ENVI_point(shapefile):
    """
    Read spatial point feature and store as a list

    :param shapefile: point feature (shapefile) with fields， such as point_x, point_y and rastervalu

    :return:
    point_list: nested list, eg, [[ID_P, X, Y, ENVI], ...[], ,,,]
    """
    point_list, rows, fields, ID = [], arcpy.SearchCursor(shapefile), arcpy.ListFields(shapefile), 10000000
    spatial_ref = arcpy.Describe(shapefile).spatialReference
    ID_P, X, Y, H = 0, 0, 0, 0
    for row in rows:
        for field in fields:
            if field.name == "POINT_X" or field.name == "point_x":
                X = row.getValue(field.name)
            if field.name == "POINT_Y" or field.name == "point_y":
                Y = row.getValue(field.name)
            if field.name == "RASTERVALU" or field.name == "rastervalu" or field.name == "GRID_CODE" or field.name == "grid_code":
                H = row.getValue(field.name)
                ID_P = "E%d" % ID  # 编号，起始点，E10000000
        point_list.append([ID_P, X, Y, H])
        ID += 1

    if len(point_list) < 1:
        raise Exception("EMPTY LIST: YOU GOT AN EMPTY LIST, PLEASE CHECK YOUR INPUT FILE!!!")
    return point_list, spatial_ref


def read_file(file_path, n=1):
    """
    This function is used to read input data from a txt file. 数据表需不包含表头，第一列中的各元素作为字典的键，同一行中余下作为值，
    其中n表示一行中，从第n各索引开始，余下部分作字典的值。

    :param file_path: 待读取文件路径
    :param n: 索引号，该位置及之后索引作为字典的值。

    :return:
    readout_dict: 数据的第一列为键，
    """
    # read file
    with open(file_path, 'r') as f:
        readout = [e.strip('\n').split(',') for e in f.readlines()]
        readout_dict = dict(zip([j[0] for j in readout], [i[n:] for i in readout]))

    # convert data type
    for k, v in readout_dict.items():
        readout_dict[k] = [float(d) for d in v]
    return readout_dict


def read_file_list(file_path, n=2):
    """
    将数据读取为列表，并将其中的数值转化为浮点数据。

    :param file_path: 待读取文件路径
    :param n: 索引号，该位置及之后需转化为浮点。

    :return:
    readout_list: 每行数据为一个子列表。
    """
    with open(file_path, 'r') as f:
        readout = [e.strip('\n').split(',') for e in f.readlines()]

    readout_list = []
    for element in readout:
        if element[0] == "0":
            continue
        else:
            x = [float(e) for e in element[n:]]
            y = element[:n]
            z = y + x
            readout_list.append(z)
    return readout_list


# ----------------------------------------------------------------------------------------------------------------------
# data pretreatment
# ----------------------------------------------------------------------------------------------------------------------


def assign_height(point_list, height, raster_size=30):
    """
    用于获取污水处理设施所处集水区的水功能区等级。

    :param raster_size: cell size of envi data set, default value is 30 meter.
    :param point_list: 污水处理设施空间坐标数据列表嵌套列表，数据结构如下，envi, level, tech, c, o&m需更新
           eg [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [envi, level],tech, c, o&m, []], ..., [], ...]
    :param height: 数值高程数据，nested list[[id, x, y, h, ...], [...], ...]

    :return:
    point_list：更新[[ID, X, Y, H, Q, M, 0, [id1, id2, id3], [envi, level],tech, c, o&m, []], ..., [], ...]中的
                envi, level, tech, c, o&m等
    """
    for i in point_list:
        X, Y = i[1], i[2]
        height_value = None
        for j in height:
            if (X + raster_size >= j[1] >= X - raster_size) and (Y + raster_size >= j[2] >= Y - raster_size):
                height_value = j[3]
                break
        i[3] = height_value
    return point_list


def check_ID(id_list, marker='H'):
    """
    用于筛选id list中的农村住户ID（前缀为“H”），其余id如DEM，ENVI，CLUSTER将被舍去。并去除列表中的重复

    :param id_list: 列表，每个元素均为point feature的ID号，不同类型point feature的ID号前缀不同。
    :param marker: ID号的前缀标识

    :return:
    out_list: 按特定前缀筛选后的id list。
    """
    out = []
    for ID in id_list:
        if ID[0] == marker:
            out.append(ID)
        else:
            continue
    out_list = list(set(out))
    out_list.sort()
    return out_list


def unfold_nested_list(nested_list):
    """
    将嵌套列表展开一级，如[[e1, e3, e3], [], ...]展开为：[e1, e2, e3, ...]
    """
    out = []
    for i in nested_list:
        out.extend(i)
    return out


def coverage_rate(point_list, rate):
    """
    根据用户输入的目标区域污水处理覆盖率，对目标区域内的households进行筛选。其筛选原则如下：删除目标区域内聚集规模较小的households，如先删除
    scatter类型的households点，再删除cluster规模为2的households，以此类推，直至满足输入的rate值。

    :param point_list: 基于输入households数据读取的数据点列表。
    :param rate: 农村污水处理覆盖率，float（%）。

    :return:
    out_list: 数据结构同point_list，只是数据点的数据有所减少
    """
    out_list = []
    point_num = len(point_list)
    point_dict = get_cluster(point_list)
    target_num = int(point_num * rate)
    delete_num = point_num - target_num  # 需删除的households数量。

    # scatter rate
    scatter_num = len(point_dict['L0'])
    cluster2_num = 0
    cluster3_num = 0
    cluster4_num = 0
    cluster5_num = 0
    for i in point_dict.values():
        if len(i) == 2:
            c2 = len(i) * 2
            cluster2_num += c2
        elif len(i) == 3:
            c3 = len(i) * 3
            cluster3_num += c3
        elif len(i) == 4:
            c4 = len(i) * 4
            cluster4_num += c4
        elif len(i) == 5:
            c5 = len(i) * 5
            cluster5_num += c5
        else:
            continue

    if delete_num == 0:  # 1
        return point_list
    elif 0 < delete_num <= scatter_num:  # 0.9
        for k, v in point_dict.items():
            if k != 'L0':
                out_list.extend(v)
            elif k == 'L0':
                rest_point = v[delete_num:]
                out_list.extend(rest_point)
    elif scatter_num < delete_num <= scatter_num + cluster2_num:  #
        list2 = []
        num = int(math.floor((delete_num - scatter_num) / 2))
        for k, v in point_dict.items():
            if k == 'L0':
                continue
            elif k != 'L0':
                if len(v) == 2:
                    list2.append(v)
                elif len(v) > 2:
                    out_list.extend(v)
        rest_point = list2[num:]
        for i in rest_point:
            out_list.extend(i)
    elif scatter_num + cluster2_num < delete_num <= scatter_num + cluster2_num + cluster3_num:
        list3 = []
        num = int(math.floor((delete_num - scatter_num) / 3))
        for k, v in point_dict.items():
            if len(v) <= 2:
                continue
            elif k != 'L0':
                if len(v) == 3:
                    list3.append(v)
                elif len(v) > 3:
                    out_list.extend(v)
        rest_point = list3[num:]
        for i in rest_point:
            out_list.extend(i)
    elif (scatter_num + cluster2_num + cluster3_num) < delete_num <= (scatter_num + cluster2_num +
                                                                      cluster3_num + cluster4_num):
        list4 = []
        num = int(math.floor((delete_num - scatter_num) / 4))
        for k, v in point_dict.items():
            if len(v) <= 3:
                continue
            elif k != 'L0':
                if len(v) == 4:
                    list4.append(v)
                elif len(v) > 4:
                    out_list.extend(v)
        rest_point = list4[num:]
        for i in rest_point:
            out_list.extend(i)
    elif (scatter_num + cluster2_num + cluster3_num + cluster4_num) < delete_num <= (scatter_num + cluster2_num + cluster3_num
                                                                                     + cluster4_num + cluster5_num):
        list5 = []
        num = int(math.floor((delete_num - scatter_num) / 5))
        for k, v in point_dict.items():
            if len(v) <= 4:
                continue
            elif k != 'L0':
                if len(v) == 5:
                    list5.append(v)
                elif len(v) > 5:
                    out_list.extend(v)
        rest_point = list5[num:]
        for i in rest_point:
            out_list.extend(i)
    return out_list


def list_to_dict(nested_list, n=0):
    """
    This function used to convert list to dictionary, ID (i[0], i is the element of list) as key, and element as value.
    This function provides a parameter (n=0) that allows the user to select the elements in the list as the key of the
    dictionary on demand. And user need to check if the chosen element meet the requirement of a dictionary key.

    :param n: the index number of element n in nested list of input list.
    :param nested_list: a nested list, eg [[id1, ...], []...]

    :return:
    dictionary: {id1: [id1, ...], ...}
    """
    key_list = [i[n] for i in nested_list]  # 嵌套列表的任意元素作为键。需留意是否重复...
    dictionary = dict(zip(key_list, nested_list))
    return dictionary


def pop_to_flow(Q, coefficient=0.08):
    """
    This function is used to convert population to sewage flow (m3).

    :param coefficient: daily sewage discharge coefficient, m3 per capita, default value is 0.08 m3/day.
    :param Q: population of considered point.

    :return:
    flow: Q * coefficient: element Q (population) is converted sewage flow (m3/d)
    """
    return Q * coefficient


def get_cluster(point_list):
    """
    This function used to convert point list to a dictionary. Marker is used as the key, the ID of point which have the
    same marker is stored in a list and make is as the value. Otherwise, the scatter points both marted as 'G0', 'L0' or
    'C0' ... so the scatter points both stored in one key-value. And the value of the dictionary is a list with details
    information of each point within each cluster, rather than id. 完整信息

    :param point_list: nested list, eg [[ID, X, Y, H, Q, M, T, [ID_LIST], [e, L], ...], ..., [], ...]

    :return:
    cluster_dict: dictionary, eg {M1,[id1, id2, ...], ..., M0: [idj,...idk]}. M0 is the scatter points.
    """
    marker = list(set([i[5] for i in point_list]))
    cluster_dict = {}
    for m in marker:
        id_list = []
        for i in point_list:
            if i[5] == m:
                id_list.append(i)
        cluster_dict[m] = id_list
        long_list = [D[0] for D in id_list]
        if len(long_list) != len(set(long_list)):
            raise Exception("len(id_list) != len(set(id_list))!!!")
    return cluster_dict


def to_one_point(near_point_list, point_dict_o):
    """
    This function is used to merge multi-point to one point. 输入的需合并的节点列表，只存储各节点的ID。

    :param near_point_list: output of get_merge_point(). id of point which needed to merge. [[node1, node2,...],  ...]
                            仅仅包含每个点的ID
    :param point_dict_o: eg {ID:[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [e, l], tech, c, o&m, [0, 0, 0, 0, 0]], ...}

    :return:
    point_list: nested point, eg [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [e, l], tech, c, o&m, []], ..., [], ...]
    """
    point_dict = deepcopy(point_dict_o)
    point_list = []
    ID = 10000
    for element in near_point_list:
        ID += 1
        ID_m = "M%d" % ID
        points = [point_dict[i] for i in element]
        x = np.mean([i[1] for i in points])
        y = np.mean([i[2] for i in points])
        h = np.mean([i[3] for i in points])
        Q = np.mean([i[4] for i in points])
        M = point_dict[element[0]][5]
        point_list.append([ID_m, x, y, h, Q, M, 0, element[:], ['e', 'l'], None, 0, 0, [0, 0, 0, 0, 0]])
    return point_list


def get_merge_point(point_list):
    """
    This function is used to merge point according to distance. 在返回值中，每个子列表列出了需合并节点的ID号。
    尽可能的考虑函数的扁平化设计。本函数只执行一轮迭代，这一轮迭代需对输入的聚集区节点进行遍历处理。

    :param point_list: 目标聚集区（散点）所包含的节点信息列表。一般而言，一个目标区域有多个目标聚集区（散点）。
                       [[ID, X, Y, H, Q, M, 0,[id1, ...], [e, l], tech, c, o&m, [0, 0, 0, 0, 0]], [], ...]

    :return:
    merged_list: nested list, 每个子列表包含带合并节点的ID号，ID为household类节点ID号。
    """
    handled_point_list = []  # id
    merged_list = []  # id
    for p in point_list:
        if p[0] not in handled_point_list:
            handled_point_list.append(p[0])
            distance_list = []  # 嵌套列表 [[distance, id], [], ...]
            for q in point_list:
                if q[0] not in handled_point_list:
                    distance = math.hypot(p[1] - q[1], p[2] - q[2])
                    distance_list.append([distance, q[0], q[7]])  # q[0]当前点p最近点的ID， q[7]该点由哪些household类节点合并而来
                else:
                    continue

            if len(distance_list) == 0:  # point_list 仅有一个点。
                merged_list.append(p[7])
            elif len(distance_list) >= 1:  # point_list有多个点（>=2）
                distance_list.sort()
                nearest_node = distance_list[0]
                handled_point_list.append(nearest_node[1])
                id_list = p[7] + nearest_node[2]
                merged_list.append(id_list)
        else:
            continue
    if len(point_list) == 0:  #
        raise Exception("EMPTY INPUTS: PLEASE CHECK YOUR INPUTS!!!")
    elif len(merged_list) == 0:
        raise Exception("EMPTY LIST: YOU GOT AN EMPTY LIST! PLEASE CHECK YOUR INPUTS!!!")
    return merged_list


def degree_of_dec(household, merged_point, charge, pop):
    """
    This function is used to calculate the degree of decentralization for rural wastewater treatment in a area.

    :param pop: 户均人口
    :param charge: 人均污水排放量
    :param household: number of households within considered area
    :param merged_point: nested point, eg [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], ...], ..., [], ...]. This list should
           contain both cluster and scatter point. 此输入数据，实际为污水处理站点的详细信息，其每个嵌套子列表的索引号为7的元素存储
           该污水处理站点所服务的农户ID号（农村居民住房空间位置数据）。

    :return:
    DD1, DD2: degree of decentralization for RuST.
    """
    W = household * charge * pop
    E = sum([i[4] * charge * len(i[7]) * len(i[7]) for i in merged_point])  # check
    DD1 = E / (W * household)
    DD2 = len(merged_point) / float(household)
    return DD1, DD2


# ----------------------------------------------------------------------------------------------------------------------
# selection of sewage treatment technology
# ----------------------------------------------------------------------------------------------------------------------


def assign_attribute(point_list, envi, raster_size=30):
    """
    用于获取污水处理设施所处集水区的水功能区等级。

    :param raster_size: cell size of envi data set, default value is 30 meter.
    :param point_list: 污水处理设施空间坐标数据列表嵌套列表，数据结构如下，envi, level, tech, c, o&m需更新
           eg [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [envi, level],tech, c, o&m, []], ..., [], ...]
    :param envi: 水功能区划数据，nested list[[id, x, y, envi], [...], ...]

    :return:
    point_list：更新[[ID, X, Y, H, Q, M, 0, [id1, id2, id3], [envi, level],tech, c, o&m, []], ..., [], ...]中的
                envi, level, tech, c, o&m等
    """
    for i in point_list:
        X, Y = i[1], i[2]
        envi_value = None
        for j in envi:
            if (X + raster_size >= j[1] >= X - raster_size) and (Y + raster_size >= j[2] >= Y - raster_size):
                envi_value = j[3]
                break
            else:
                envi_value = 5
        i[8][0] = envi_value
    return point_list


def get_standard(facility_list, standard, coefficient=0.08):
    """
    根据污水处理设施规模，集水区水功能区等级确定所需达到的污水排放标准。

    :param facility_list: facility list with field  ENVI. output of function assign_attribute.
           eg: [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [envi, level],tech, c, o&m, []], ..., [], ...],
           其中level, tech, c, o&m需更新, Q任然为户均人口
    :param coefficient: daily sewage discharge coefficient, m3 per capita, default value is 0.08 m3/day.
    :param standard: 地方/国家农村WWTP出水排放标准

    :return:
    facility_list: updated list, field envi was added.
    """
    for e in facility_list:  # todo check
        scale = math.ceil(len(e[7]) * pop_to_flow(e[4], coefficient))  # covert population to sewage flow, 80L/p.day
        envi = int(e[8][0])
        for _, s in standard.items():
            if s[0] <= scale <= s[1] and envi == s[2]:
                e[8][1] = int(s[3])
            else:
                continue
    return facility_list


def system_option(facility_list, technology_list, charge=0.08, facility_ls=30):
    """
    used to select the most cost-effective wastewater treatment technology for each facility within research area.

    :param facility_list: facility list with field  ENVI. output of function assign_attribute.
           eg: [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [envi, level], tech, c, o&m, []], ..., [], ...],
           其中tech, c, o&m需更新, Q任然为户均人口
    :param technology_list: with fields ['id', 'name', 'level', 'scale_min', 'scale_max', 'cost', 'o&m']
    :param charge: daily sewage discharge coefficient, m3 per capita, default value is 0.08 m3/day.
    :param facility_ls: lifespan of sewage treatment facility

    :return:
    facility_list: facility_list with technology and price (construction, operation and maintenance)
                   eg. [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [envi, level], tech, c, o&m], ..., [], ...]
    """

    def sort_index(elem):
        """
        按嵌套列表中每个嵌套元素索引号为5的元素排序。该函数作为sort（)函数中的key值使用。

        :param elem: ['id', 'name', 'level', 'scale_min', 'scale_max', 'cost', 'o&m']
        """
        return elem[5]

    for e in facility_list:
        scale = math.ceil(len(e[7]) * pop_to_flow(e[4], charge))  # covert population to sewage flow, 80L/p.day
        standard = e[8][1]
        tec_list = []
        # find suitable technology
        for i in technology_list:
            price, level = i[5], i[2]
            scale_min, scale_max = i[3], i[4]
            if level <= standard and scale_min <= scale <= scale_max:
                tec_list.append(i)

        # find the most cost-effective technology
        if len(tec_list) == 0:
            raise Exception("ERROR: NO SUITABLE TECHNOLOGY!!!")
        else:
            tec_list.sort(key=sort_index)
            min_cost = tec_list[0]
            e[9], e[10], e[11] = min_cost[1], min_cost[5], min_cost[6] * scale * facility_ls  # 更新嵌套列表（造价、运维费用）
    return facility_list


def pipe_cost_onsite(point_list, price_list, lifespan=30):
    """
    This function is used to calculate the cost of sewage collection system.
    用于计算每套污水处理设施所配套管网系统的造价、运维费用。最终返回目标区域所有系统配套管网设施的造价、运维费用。管网系统的设计寿命为30年

    :param lifespan: 使用寿命
    :param point_list: nested point, eg [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], ...], ..., [], ...]，该列表还包含了散户。
    :param price_list: 不同规模[p1, p2, p3, p4, p5, p6]设施配套管网系统的户均造价，包括了管道、检查井、化粪池、隔油井，泵站等。

    :return:
    sewer_c, sewer_m: 化粪池和隔油井造价，管网系统的建设费用、运维费用。本函数输出为目标区域内的总造价
    """
    sewer_c, sewer_m = 0, 0  # 化粪池和隔油池的综合造价定为3500元
    for n in point_list:
        num = len(n[7])
        for p in price_list:
            floor, upper, price = p
            if floor <= num <= upper:
                cost = num * price
                sewer_c += cost
            else:
                continue
    sewer_m = lifespan * sewer_c * 0.011
    return sewer_c, sewer_m


# ----------------------------------------------------------------------------------------------------------------------
# Merging model
# ----------------------------------------------------------------------------------------------------------------------


def get_distance(source, sink, dis=50):
    """
    This function is used to calculate the three-dimensional Euclidean distance of two points. And check how many
    inspection well is needed.

    :param dis: distance between two inspection well, default value is 50 meter.
    :param source: start point, list [ID, X, Y, H, Q, M, 0.5, [ID], [envi, level], tech, c, o&m, []].
    :param sink: end point, list [ID, X, Y, H, Q, M, 0.5, [ID], [envi, level], tech, c, o&m, []].

    :return:
    distance: Euclidean distance of the given two points, float
    distance_2D: 2D Euclidean distance of he given two points, float
    inspection: number of inspection well, int.
    """
    height = abs(source[3] - sink[3])
    distance_2D = np.hypot((source[1] - sink[1]), (source[2] - sink[2]))
    distance = np.hypot(height, distance_2D)
    inspection = int(distance_2D / dis) + 1
    return distance, distance_2D, inspection


def get_diameter(Q, type_list):
    """
    This function is used to ge the diameter of pipeline. And converted to standard specifications. 本项目中使用管径规格如
    下：

    :param Q: Flow of sewage, float (m3).
    :param type_list: dictionary (ascending order) with standard specifications of pipe, (mm)

    :return:
    diameter: standard specifications of pipe
    """
    diameter = None
    for i in sorted(type_list.keys()):
        Q_max = 0.8 * 9 * 3600 * 3.14159 * (int(i) ** 2) / 4000000  # m3/day, flow speed is 1 m/s
        if Q <= Q_max:
            diameter = int(i)
            break
        else:
            continue
    if None == diameter:  # raise error, if the diameter in type_list do not meet the requirement of Q.
        raise Exception("Please input more larger pipe specification!!!")
    return diameter


def get_earthwork(source, sink, diameter, depth_min=0.5, wide=0.5, depth_max=4, slope_d=0.005):
    """
    This function is used to calculate the volume of earthwork. If the slope of the pipeline is greater than the default
    value a drop well is required. The flow direction is from source point to sink point.
    输入数据中的source和sink点可能是数据结构相同的各类point feature，模型初始读取的point feature数据中，默认索引号为7的子元素（子列表）
    存储改点自身ID号，在后续计算中该位置存储ID号应为实际产生污水的的。因此，对于本身无污水产生的点（如DEM）需从数据中剔除。

    :param wide: wide of trench, meter
    :param diameter: diameter of pipeline, mm
    :param source: start point of a pipeline; [ID, X, Y, H, Q, M, 0.5,[id1, id2, id3], [envi, level], tech, c, o&m, []].
    :param sink: end point of a pipeline; [ID, X, Y, H, Q, M, 0.5,[id1, id2, id3], [envi, level], tech, c, o&m, []].
    :param depth_min: the minimum depth of pipeline, the default value is 0.5 meter.
    :param depth_max: the maximum depth of pipeline, the default value is 4 meter.
    :param slope_d: the slope of pipeline and the default value is 0.005. 管道设计坡度

    :return:
    earthwork: the volume of earthwork, float type(m3).
    drop: 0 or 1. 1 means a drop well is required; 0 means no drop well.
    pump_station: number of pump station
    """
    earthwork, drop, inspection, pump_station = None, 0, 0, 0
    w = diameter / 1000.0 + wide  # meter
    height = source[3] - sink[3]  # 高程差
    trench = sink[6]  # 终点当前埋深，如果该点有另外支路接入，需终点考虑。
    distance, distance_2d, _ = get_distance(source, sink)
    depth = depth_min + distance_2d * slope_d - height  # 埋深，理论计算值。
    # calculate earthwork
    if height >= 0:  # 顺坡：①小于最小值，②范围内，③大于最大值
        if depth > depth_max:
            pump_station = 1
            inspection = 1
            earthwork = 2 * 0.5 * depth_min * distance_2d * w  # 开挖和回填
            household_nodes = sink[7] + source[7]  # flow
            pop = check_ID(household_nodes, "H")
            sink[7] = pop
            sink[6] = depth_min  # trench
        elif depth <= depth_max:  # trench = 0.5-4
            if depth > trench:
                earthwork = 2 * 0.5 * (source[6] + depth) * distance_2d * w
                inspection = 1
                household_nodes = sink[7] + source[7]  # flow
                pop = check_ID(household_nodes, "H")
                sink[7] = pop
                sink[6] = depth  # trench
            elif depth_min <= depth <= trench:
                earthwork = 2 * 0.5 * (source[6] + depth) * distance_2d * w
                inspection = 1
                drop = 1
                household_nodes = sink[7] + source[7]  # flow
                pop = check_ID(household_nodes, "H")
                sink[7] = pop
                sink[6] = trench  # trench
            elif depth < depth_min:
                earthwork = 2 * 0.5 * (source[6] + depth_min) * distance_2d * w
                drop = 1
                household_nodes = sink[7] + source[7]  # flow
                pop = check_ID(household_nodes, "H")
                sink[7] = pop
                sink[6] = depth_min  # trench
    else:  # height < 0:  逆坡：①范围内，②大于最大值
        if depth > depth_max:
            pump_station = 1
            inspection = 1
            earthwork = 2 * depth_min * distance_2d * w
            household_nodes = sink[7] + source[7]  # flow
            pop = check_ID(household_nodes, "H")
            sink[7] = pop
            sink[6] = depth_min  # trench
        elif depth <= depth_max:
            if depth > trench:
                inspection = 1
                earthwork = 2 * 0.5 * (source[6] + depth) * distance_2d * w
                household_nodes = sink[7] + source[7]  # flow
                pop = check_ID(household_nodes, "H")
                sink[7] = pop
                sink[6] = depth  # trench
            elif depth_min <= depth <= trench:
                earthwork = 2 * 0.5 * (source[6] + depth) * distance_2d * w
                inspection = 1
                drop = 1
                household_nodes = sink[7] + source[7]  # flow
                pop = check_ID(household_nodes, "H")
                sink[7] = pop
                sink[6] = trench  # trench
    return earthwork, drop, inspection, pump_station


def pipeline_cost(earthwork, pipeline, diameter, drop, inspection, pump, price_list, lifespan=30):
    """
    This function is used to calculate the total cost of pipeline, contains construction cost( material and earthwork),
    operation and maintenance cost. 两节点间的管道造价，配合get_earthwork使用。
    不包括化粪池、隔油井的相关费用。

    :param lifespan: life span of pipeline, default value is 30 year.
    :param pump: number of pump station, int
    :param earthwork: excavation and back-filling, m3
    :param pipeline: 3D Euclidean length of pipeline, meter.
    :param diameter: pipeline diameter, mm.
    :param drop: 0 or 1, 1 means a drop well is required; 0 means no drop well.
    :param inspection: number of inspection well, int
    :param price_list: price of earthwork, pipe, drop well, inspection well, pump,eg[12,{type: 59, ...}, 32, 44, 2344]


    :return:
    cost: total cost of pipeline, cny
    cost_con: construction cost, cny.
    cost_om: operation and maintenance (o&m) cost, cny.
    """
    pipe_price_dict = price_list[1]  # type is set as key, eg. {100: 34}
    cost_con = earthwork * price_list[0][0] + pipeline * pipe_price_dict[str(diameter)][0] + \
               drop * price_list[2][0] + inspection * price_list[3][0] + pump * price_list[4][0]
    cost_om = 0.01 * cost_con * lifespan
    cost = cost_con + cost_om
    return cost, cost_con, cost_om


def get_h_value(point, end, price=150, lifespan=30):
    """
    This function is used to estimate the cost from point x to the end point, use the Manhattan distance.
    A-star 启发式函数H值的计算，用于估计中间点至终点的代价。

    :param point: the coordinate of point x; eg [ID, x, y, z, q, ...].
    :param end: the coordinate of end point; eg [ID, x, y, z, q, ...].
    :param price: comprehensive cost of pipeline, default value is 150 CNY。不含检查井
    :param lifespan: life span of pipeline, default value is 30 year.

    :return:
    h_value: 中间点x至终点（目标点）的估计代价。
    """
    distance = abs(point[1] - end[1]) + abs(point[2] - end[2]) + abs(point[3] - end[3])
    con = distance * price + (1 + int(distance / 50)) * 3000
    h_value = con + 0.01 * con * lifespan
    return h_value


def get_neighbor(point, graph, coefficient=1.5, cellsize=30):
    """
    This function is used to get the neighbor points within the GRAPH of a given point P.
    用于获取底图中与目标点（当前点）相邻的各点（包括斜对角点）

    :param point: the coordinate of point x; [ID, X, Y, H, Q, M, 0.5,[id1, id2, id3], [envi, level], tech, c, o&m, []].
    :param graph: point list, points with coordinate; [ID, X, Y, H, Q, M, 0.5,[id], [envi, level], tech, c, o&m, []]
                  graph实际就是DEM point feature数据读取后的point list。
    :param cellsize: raster size, 30 meter is set as default value
    :param coefficient: used to calculate the scanning distance of the given point P, default value is 1.5.

    :return:
    neighbor_dict: dictionary with neighbor points of point P.
    eg: {ID: [[ID, X, Y, H, Q, M, 0.5,[id], [envi, level], tech, c, o&m, [ID, farther, F, G, H]], ...}
    """
    neighbor_dict = {}
    near_distance = cellsize * coefficient
    for i in graph:
        distance = np.hypot(point[1] - i[1], point[2] - i[2])
        if distance <= near_distance and distance != 0:
            neighbor_dict[i[0]] = i + [[i[0], point[0], 0, 0, 0]]  # todo None or 0

    if len(neighbor_dict) == 0:
        raise ValueError("You got AN EMPTY neighbor list! Please check you inputs!!!")
    return neighbor_dict


def get_weight(start, point, end, diameter, price_list, depth_min=0.5, wide=0.5, depth_max=4, slope_d=0.005,
               price=150, lifespan=30, dis=50):
    """
    This function is used to calculate or upgrade F, G and H value. 计算中间点的启发式函数权重。这个中间点为当前处理点的相邻点。

    :param price_list: price of earthwork, pipe, drop well, inspection well, pump,eg[12,{type: 59, ...}, 32, 44, 2344]
    :param dis: distance between two inspection well, default value is 50 meter.
    :param depth_min: the minimum depth of pipeline, the default value is 0.5 meter.
    :param depth_max: the maximum depth of pipeline, the default value is 4 meter.
    :param slope_d: the slope of pipeline and the default value is 0.005. 管道设计坡度
    :param wide: wide of trench, meter
    :param diameter: diameter of pipeline, mm
    :param start: the front point of current, list [id, X, Y, H, Q, M, t,[id], [], tech, c, o&m, []].
    :param point: point need to calculate F, G and H. [id, X,Y,H,Q, M, t,[id], [], tech, c, o&m, [ID, farther, F, G, H]]
    :param end: the back point of current, list [id, X, Y, H, Q, M, t,[id], [], tech, c, o&m, []].
    :param price: comprehensive cost of pipeline
    :param lifespan: ife span of pipeline, default value is 30 year.

    :return:
    point: [id, X,Y,H,Q, M, t,[id], [], tech, c, o&m, [ID, farther, F, G, H]] updated
    """
    pipeline, _, _ = get_distance(start, point, dis)
    earthwork, drop, inspection, pump = get_earthwork(start, point, diameter, depth_min, wide, depth_max, slope_d)
    G_value, _, _ = pipeline_cost(earthwork, pipeline, diameter, drop, inspection, pump, price_list, lifespan)
    H_value = get_h_value(point, end, price, lifespan)
    point[13][2], point[13][3], point[13][4] = G_value + H_value, G_value, H_value
    return point


def get_path(start, end, path_dic):  # 核实路径用顺序好还是逆序好 todo
    """
    This function is used to get the path from start point to end point.
    返回值记录起点至终点最优路径所经过的每一个节点的ID号，并按经过的先后逆序排列。

    :param start: ID of start point.
    :param end: ID of end point.
    :param path_dic: point ID dictionary, {ID: [ID, parent],...}。路径字典

    :return:
    path: pipeline path from source point to sink point. [start, id1, id2, ... source]
    """
    path, check = [end], end
    while check:   # or check is not None:
        for v in path_dic.values():
            if check == v[0]:
                if v[1] is None:
                    check = v[1]
                else:
                    path.append(v[1])
                    check = v[1]
            else:
                continue

    path.reverse()
    return path


def a_star(source, sink, graph, type_list, price_list, coeff=1.5, depth_min=0.5, wide=0.5, depth_max=4, slope=0.005,
           cellsize=30, lifespan=30, dis=50):
    """
    This function is used to calculated the path from start point to end point based on A* algorithm. The weight used to
    find the optimal path (most cost-effective) is the total cost of pipeline, rather the distance.

    Do not consider the loss of flow.

    :param lifespan: ife span of pipeline, default value is 30 year.
    :param wide: wide of trench, meter
    :param price_list: price of earthwork, pipe, drop well, inspection well, pump,eg[12,{type: 59, ...}, 32, 44, 2344]
    :param source: point coordinate of start point, [id, X,Y,H,Q, M, t,[id], [], tech, c, o&m, [ID, '', 0, 0, 0]].
    :param sink: point coordinate of end point, [ID, x, y, z, q].
    :param graph: point list，eg [[id, X,Y,H,Q, M, t,[id], [], tech, c, o&m, [ID, farther, F, G, H]], [], ...]
    :param depth_min: the minimum depth of pipeline, default value is 0.5 meter.
    :param depth_max: the maximum depth of pipeline, default value is 4 meter.
    :param slope: the slope of pipeline, default value is 0.005.
    :param coeff: used to calculate the scanning distance of the given point P, default value is 1.5.
    :param type_list: list (ascending order) with standard specifications of pipe, (mm)
    :param cellsize: raster size, 30 meter is set as default value
    :param dis: distance between two inspection well, default value is 50 meter.

    :return:
    path: [(start, end), cost, [ids, id1, id2, ..., idn, ide]]
    node_dict: details for every node of path.
               {(ids, ide): [[id, X,Y,H,Q, M, t,[id], [], tech, c, o&m, [ID, '', 0, 0, 0]], [], ...]}
    """
    # data pretreatment
    H = get_h_value(source, sink)
    path, cost, close_dict = [], None, {}
    open_list = [[0, 0, H, source[0], None]]  # F, G, H, POINT, PARENT
    flow_to_sink = source[4] * len(source[7]) * 0.08
    diameter = get_diameter(flow_to_sink, type_list)
    graph_node = deepcopy(graph)  # 引用传递，导致数据污染
    graph_node.append(source)
    graph_node.append(sink)
    graph_dict = list_to_dict(graph_node)
    # get neighbor points
    while open_list:
        open_list.sort()  # 确保在原位上修改而非产生一个新对象（sorted())
        parent = open_list.pop(0)  # [F, G, H, POINT, PARENT]
        if parent[3] == sink[0]:  # it is necessary to ensure the ID number keep constant
            close_dict[parent[3]] = parent
            cost = parent[0]
            break
        point = graph_dict[parent[3]]  # parent point
        close_dict[parent[3]] = parent
        neighbors = get_neighbor(point, graph_node, coeff, cellsize)
        for k, P in neighbors.items():  # get F, H, G (当前点、父节点、终点，而非当前点、起点、终点）
            if P[0] in close_dict:
                continue
            else:
                pnt = get_weight(point, P, sink, diameter, price_list, depth_min, wide,
                                 depth_max, slope, 70, lifespan, dis)
                PNT = pnt[13]  # get nested list [id, father, F, G, H]
                P = [PNT[0], PNT[1], parent[1] + PNT[3] + PNT[4], parent[1] + PNT[3], PNT[4]]
                if k in [i[3] for i in open_list]:  # update weight
                    index_num = None
                    for i in open_list:
                        if k == i[3]:
                            index_num = open_list.index(i)
                    if P[3] < open_list[index_num][1]:  # more cost-effective from current point to the neighbor
                        open_list[index_num] = [P[2], P[3], P[4], P[0], P[1]]
                    else:  # no need to update
                        continue
                else:
                    open_list.append([P[2], P[3], P[4], P[0], P[1]])
    # get path
    point_dict = {}
    for k, v in close_dict.items():
        point_dict[k] = v[3:]
    path_list = get_path(source[0], sink[0], point_dict)
    path = [(source[0], sink[0]), cost, path_list]
    return path


# ----------------------------------------------------------------------------------------------------------------------
# A greedy algorithm used for optimal sewage pipeline designing.
# ----------------------------------------------------------------------------------------------------------------------


def get_scenario_cluster(scenario_list, marker='L0'):
    """
    以simulate_to_centralized函数输入scenario_list为输入，获取不同情景下每个WWTP所服务的cluster或scatter类点

    :param scenario_list: cluster point list, communities or scatters need to be merged is stored in a list and element
                          within this list represents a node. cluster_list [[node1, node2, ...], [], ...]. and the data
                          structure of a node is as follow:
                          [id, X, Y, H, Q, M, t,[id, ...], [], tech, c, o&m, [ID, farther, F, G, H]]
    :param marker: households 节点中散点的聚类标识。 默认为" L0"

    :return:
    scenario_cluster: nested list, [[cluster1, scatter1, ...], [], ...]. and the data structure of a node is as follow:
                      [id, X, Y, H, Q, M, t,[id, ...], [], tech, c, o&m, [ID, farther, F, G, H]]
    nodes_list: 嵌套列表，存储目标区域内所有cluster or scatter类节点。
    """
    def s(e):
        return e[7]

    scenario_cluster, nodes_list = [], []

    for cluster in scenario_list:  # mubiao quyu de yige wushui chulizhandian de fuwu duixiang (scatter or cluster).
        marker_list = list(set([i[5] for i in cluster]))
        cluster_list = []
        cluster_dict = {}
        for m in marker_list:
            node = []
            for h in cluster:  # 每个污水处理站点所服务的households
                if h[5] == m:
                    node.append(h)
                else:
                    continue
            cluster_dict[m] = node

        for k, v in cluster_dict.items():
            if k == marker:
                cluster_list.extend(cluster_dict[k])
            else:  # k != 'L0'
                X = np.mean([x[1] for x in v])
                Y = np.mean([y[2] for y in v])
                H = np.mean([h[3] for h in v])
                Q = np.mean([q[4] for q in v])
                M = v[0][5]
                id_list = [i[7] for i in v]
                ID_list = []
                for i in id_list:
                    ID_list.extend(i)
                    ID_list.sort()
                cluster_list.append([0, X, Y, H, Q, M, 0, list(set(ID_list)), ['e', 'l'], 't', 0, 0, [0, 0, 0, 0, 0]])

        scenario_cluster.append(cluster_list)
        nodes_list.extend(cluster_list)

    # get ID
    cnt = 500000
    key_list = []
    nodes_list.sort(key=s)
    for i in nodes_list:
        if i[5] == marker:
            continue
        else:
            cnt += 1
            ID = "M%d" % cnt
            i[0] = ID
            k = sorted(i[7])
            key_list.append(tuple(k))
    nodes_dict = dict(zip(key_list, nodes_list))

    # 更新scenario_cluster中的ID号
    for C in scenario_list:
        for K, V in nodes_dict.items():
            for H in C:
                if list(K) == H[7].sort():
                    H[0] = V[0]  # ID
                else:
                    continue
    return scenario_cluster, nodes_list


def get_sewer_system(scenario_list, graph, type_list, price_list, coeff=1.5, depth_min=0.5, wide=0.5,
                     depth_max=4, slope=0.005, cellsize=30, lifespan=30, dis=50):
    """
    This function is used to generate the most cost-effective sewage collection system of several rural residential
    points. We assume that WWTP is built in the node which with lowest elevation. The final results of this function is
    a dictionary with the sewer network.

    :param scenario_list: cluster point list, communities or scatters need to be merged is stored in a list and element
                          within this list represents a node. cluster_list [[node1, node2, ...], [], ...]. and the data
                          structure of a node is as follow:
                          [id, X, Y, H, Q, M, t,[id, ...], [], tech, c, o&m, [ID, farther, F, G, H]]
    :param graph: nested list [[id, X, Y, H, 0, o, t,[id, ...], [], tech, c, o&m, [ID, farther, F, G, H]], [], ... ]
    :param type_list: type_list: list (ascending order) with standard specifications of pipe, (mm)
    :param price_list:  price of earthwork, pipe, drop well, inspection well, pump,eg[12,{type: 59, ...}, 32, 44, 2344]
    :param coeff: used to calculate the scanning distance of the given point P, default value is 1.5.
    :param depth_min: the minimum depth of pipeline, default value is 0.5 meter.
    :param depth_max: the maximum depth of pipeline, default value is 4 meter.
    :param wide: wide of trench, the default value is 0.5 meter
    :param slope: the slope of pipeline, default value is 0.005.
    :param cellsize: raster size, 30 meter is set as default value
    :param lifespan: ife span of pipeline, default value is 30 year.
    :param dis: distance between two inspection well, default value is 50 meter.

    :return:
    sewer_c, sewer_o: 分别表示目标区域内管道系统的造价和运维费用（不包括户间管网系统）
    output_path_dict：情景x的需合并污水处理站点间的配套管网数据{(sink):[[path], [path], ... ], ...}，其中path为A-star算法输出。
                      sink包含node的详细信息，非ID号。
    """

    def sort3(ele):
        """sort list according to the value of the 4th element"""
        return ele[3]

    # get cluster
    scenario_cluster, nodes_list = get_scenario_cluster(scenario_list, marker='L0')
    scenario_cluster, nodes_list = deepcopy(scenario_cluster), deepcopy(nodes_list)  # deepcopy
    point_list = nodes_list + graph
    point_dict = list_to_dict(point_list)
    output_path_list = []
    cnt = 0
    # 处理每个污水处理站点的配套管网
    for element in scenario_cluster:  # cluster中的元素为信息完整的点，而非ID
        if len(element) == 1:
            cnt += 1
            output_path_list.append([[(element[0][0], element[0][0]), 0, []]])
        elif len(element) == 2:  # 知道两点间的管网造价即可
            cnt += 1
            element.sort(key=sort3)
            sink, source = element
            path = a_star(source, sink, graph, type_list, price_list, coeff, depth_min, wide, depth_max, slope,
                          cellsize, lifespan, dis)
            output_path_list.append([path])

        elif len(element) > 2:
            # data pretreatment. 获取WWTP位点（sink），以及sink点以外的聚集点和散点。
            element.sort(key=sort3)  # 包括完整信息
            sink = element[0]
            path_node_list = []  # 用于存储路径节点，包括核心节点
            path_list = []  # 存储两点间的管网造价、路径，a*输出。
            handled_n = [sink]  # 已经处理后的点，包括sink   # 完整信息，非ID
            while len(element) > len(handled_n):
                cnt += 1
                if len(path_node_list) == 0:  # 初始运算时(最远点），算法与后续迭代有所不同（最近点）
                    cluster_node_dict, path_node_dict = {}, {}
                    for E in element:  # find the node with the longest distance to sink. 所有非sink之外的node
                        if E not in handled_n:  # 找出距离最远的点（曼哈顿距离）
                            distance_m = abs(sink[1] - E[1]) + abs(sink[2] - E[2]) + abs(sink[3] - E[3])
                            if distance_m not in cluster_node_dict.keys():
                                cluster_node_dict[distance_m] = E
                            else:
                                continue
                        else:
                            continue

                    distance = sorted(cluster_node_dict.keys(), reverse=True)
                    source_node = cluster_node_dict[distance[0]]
                    handled_n.append(source_node)
                    path = a_star(source_node, sink, graph, type_list, price_list, coeff, depth_min, wide, depth_max,
                                  slope, cellsize, lifespan, dis)
                    path_list.append(path)
                    path_node_list.extend(path[2])

                else:  # len(path_node_list) != 0
                    distance_dict = {}
                    for C in element:
                        path_dict_c = {}
                        if C not in handled_n:  # 各点至管网最近的距离
                            for node in path_node_list:  # get nearest path node to, ID
                                node_inf = point_dict[node]
                                distance_mc = abs(C[1] - node_inf[1]) + abs(C[2] - node_inf[2])
                                if distance_mc not in path_dict_c.keys():
                                    path_dict_c[distance_mc] = [C, node_inf]  # [SOURCE, SINK]
                                else:
                                    continue
                        else:
                            continue
                        distance_list = sorted(path_dict_c.keys())  # nearest distance
                        distance_dict[distance_list[0]] = path_dict_c[distance_list[0]]
                    max_distance = sorted(distance_dict.keys(), reverse=True)  # longest distance
                    source_node, sink_node = distance_dict[max_distance[0]]
                    handled_n.append(source_node)
                    path = a_star(source_node, sink_node, graph, type_list, price_list, coeff, depth_min, wide,
                                  depth_max, slope, cellsize, lifespan, dis)
                    path_list.append(path)
                    path_node_list.extend(path[2])
            output_path_list.append(path_list)
            gc.collect()
    return output_path_list


def from_node_filter(point_list):
    """
    用于检查管网节点中的上游节点是否产生污水，即区分该节点是否为household点，该类点有污水排出，而DEM点无生活污水排放。计算WWTP
    污水处理规模、管段管径/流量时，需确定汇入的household数目。过滤掉记录中不产生污水的节点ID号。

    :param point_list: nested list [[id, X, Y, H, 0, o, t,[id, ...], [], tech, c, o&m, []], [], ... ]

    :return:
    point_list: 剔除[id, ...]中不nested list [[id, X, Y, H, 0, o, t,[id, ...], [], tech, c, o&m, []], [], ... ]
    """
    for point in point_list:
        p_id_list = []
        for p_id in point[7]:
            if p_id[0] == "H":  # household 类节点ID前缀"H"
                p_id_list.append(p_id)
            else:
                continue
        point[7] = p_id_list  # update
    return point_list


def get_sewer_nodes(path_list, graph, cluster_point):
    """
    This function is used to obtain the nodes of sewer network and each node contains detail information. This function
    also return the WWTP list within the considered area.

    :param path_list: sewer path connect sink point(WWTP or node of sewer network).
    :param graph: nested list [[id, X, Y, H, 0, M, D,[id, ...], [], tech, c, o&m, [ID, farther, F, G, H]], [], ... ]
    :param cluster_point: eg [[ID, X, Y, H, Q, M, D, [id,...], [, ], tech, c, o&m], ..., [], ...], 0 为管底标高。
           and H, Q, M need to be updated, [id, ...] stores all households id within this cluster.

    :return:
    path_node_list: This is a nested list. [[id, x, y, z, Q, M, D, [id, id...]], [], ...], z is the value of elevation,
                    Q is population,  D is trench depth.
    WWTP_list: WWTP列表，含详细信息，[ENVI, level], tech, c, o&m等需更新
    """
    # pretreatment
    path_node_list = []
    point_list = graph + cluster_point  # graph中Q=‘Q’， cluster_point中Q=pop(2.7)
    point_dict = list_to_dict(point_list)
    path_node_id = []
    cluster_dict = list_to_dict(cluster_point)

    # get sewer nodes
    for network in path_list:  # check todo
        if len(network) > 0:
            for path in network:
                path_node_id.extend(path[2])
        else:
            continue
    path_node_id = list(set(path_node_id))  # 去重复
    for node_id in path_node_id:
        node = point_dict[node_id]
        path_node_list.append(node)
    path_node_list = from_node_filter(path_node_list)

    # get wwtp
    WWTP_list = []  # WWTP_LIST [[[id, X, Y, H, 0, M, T,[id, ...], [], tech, c, o&m, []], ...]
    for network in path_list:  # network为一个WWTP的配套管网系统, eg. [[('H100034', 'H100033'), 8478.9619383052959, ['H100034', 'H100033']], ...]
        for pipe in network:  # network中的一个管段
            WWTP_id = network[0][0][-1]
            WWTP = cluster_dict[WWTP_id]
            if pipe[2]:   # 非空（如[]或None）
                for n in pipe[0]:
                    if n in cluster_dict:
                        WWTP[7].extend(cluster_dict[n][7])
                    else:
                        continue
                    WWTP[7] = list(set(WWTP[7]))
                WWTP_list.append(WWTP)
            else:
                WWTP_list.append(WWTP)
    return path_node_list, WWTP_list


def get_sewer_path(path_list):
    """
    This function is used to generate the pipeline path of each sewage treatment facility

    :param path_list: sewer path connect sink point(WWTP or node of sewer network), 以wwtp点作为字典的键，{(sink): [path]}
                      record the sewer nodes' id. 目标区域内所有的wwtp

    :return:
    sewer_path_dict: sewer path of a facility(WWTP or node of sewer network), 以wwtp点作为字典的键，{(sink): [path]}
                     record the sewer pipeline(start_id, end_id)
    source_node_dict:  Source nodes of a facility. {WWTP1: [ID1, ID2, ...], ...}   ** ID
    """
    sewer_path_dict = {}
    source_node_dict = {}
    for path_l in path_list:  # pipe network system of a facility
        WWTP = path_l[0][0][1]  # ID number
        sewer_path_list = []
        source_node_list = [WWTP]  # WWTP[0]
        for pipe in path_l:  # one pipeline of pipe network system of a facility
            if pipe[2]:  # 非空
                path_node_from = pipe[2][:-1]  # 起始端
                path_node_to = pipe[2][1:]  # 终端
                sewer_path = zip(path_node_from, path_node_to)
                sewer_path_list.extend(sewer_path)
                source_node_list.append(pipe[0][0])
            else:
                source_node_list = [WWTP]
            sewer_path_dict[WWTP] = sewer_path_list  # {WWTP1: [(), (), ...], ...}
            source_node_dict[WWTP] = source_node_list
    return sewer_path_dict, source_node_dict


def get_collection_system_cost(path_list, graph, cluster_point, pipe_type, price_list, charge=0.08, depth_min=0.5,
                               wide=0.5, depth_max=4, slope=0.005, lifespan=30):
    """
    This function is used to get the design parameters (Flow/Q, trench, diameter) of pipeline. And finally return the
    total construction of sewage collection system. 每个聚集区内部管网、化粪池、隔油井造价不含在内。

    :param path_list: sewer path connect sink point(WWTP or node of sewer network), 以sink点作为字典的键。包含目标区域内所有
                      污水处理站点的的配套管网，数据以每个污水处理站点为存储单元，每个污水处理站点的管网系统包括若干（0,1,2...）管段
                      eg.[path1, path2, ...], path1: [(start, end), cost, [start, ..., end]
    :param graph: nested list [[id, X, Y, H, 0, o, t,[id, ...], [], tech, c, o&m, [ID, farther, F, G, H]], [], ... ]
    :param cluster_point: eg [[ID, X, Y, H, Q, M, 0, [id,...], [, ], tech, c, o&m], ..., [], ...], 0 为管底标高。
           and H, Q, M need to be updated, [id, ...] stores all households id within this cluster.
    :param price_list: price of earthwork, pipe, drop well, inspection well, pump,eg[12,{type: 59, ...}, 32, 44, 2344]
    :param charge: daily sewage discharge coefficient, m3 per capita, default value is 0.08 m3/day.
    :param pipe_type: list (ascending order) with standard specifications of pipe, (mm)
    :param depth_min: the minimum depth of pipeline, the default value is 0.5 meter.
    :param wide: wide of trench, meter
    :param depth_max: the maximum depth of pipeline, the default value is 4 meter.
    :param slope: the slope of pipeline and the default value is 0.005. 管道设计坡度
    :param lifespan: life span of pipeline, default value is 30 year.

    :return:
    pipeline_cost_dict: 目标区域内，cluster node间铺设管网系统的造价。每个聚集区内部管网、化粪池、隔油井造价不含在内。
    path_node_dict: path nodes with details information
    """
    # data pretreatment
    path_node_list, _ = get_sewer_nodes(path_list, graph, cluster_point)  # 包含不同情景下的所有WWTP和管网节点
    path_node_dict = list_to_dict(path_node_list)
    sewer_path_dict, source_node_dict = get_sewer_path(path_list)
    pipeline_cost_dict = {}  # 结果输出

    # iteration each WWTP
    for WWTP in sewer_path_dict:  # 遍历每个污水处理站点（WWTP）#
        path = sewer_path_dict[WWTP]  # all pipeline path of a WWTP
        source_nodes = source_node_dict[WWTP]  # ID
        sewer_cost, sewer_cost_con, sewer_cost_om = 0, 0, 0  # total cost, construction cost and O&M
        if len(path) == 0:  # NO SEWER SYSTEM
            pipeline_cost_dict[WWTP] = [sewer_cost, sewer_cost_con, sewer_cost_om]
        else:
            handled_pipe = []
            # 迭代path，直至WWTP
            while len(handled_pipe) - len(path) < 0:  # ID !=
                leaves_nodes = source_nodes[:]
                next_nodes_list = []
                for pipe in path:  # 迭代管段
                    if pipe[0] in leaves_nodes: # and pipe not in handled_pipe:
                        from_node = path_node_dict[pipe[0]]
                        to_node = path_node_dict[pipe[1]]
                        handled_pipe.append(pipe)
                        next_nodes_list.append(pipe[1])
                        distance, _, _ = get_distance(from_node, to_node)
                        flow = pop_to_flow(from_node[4] * len(from_node[7]), charge)
                        diameter = get_diameter(flow, pipe_type)
                        earthwork, drop, inspection, pump = get_earthwork(from_node, to_node, diameter, depth_min, wide,
                                                                          depth_max, slope)
                        cost, cost_con, cost_om = pipeline_cost(earthwork, distance, diameter, drop, inspection, pump,
                                                                price_list, lifespan)
                        sewer_cost += cost
                        sewer_cost_con += cost_con
                        sewer_cost_om += cost_om
                    else:
                        continue

                pipeline_cost_dict[WWTP] = [sewer_cost, sewer_cost_con, sewer_cost_om]
                source_nodes = next_nodes_list[:]  # shallow copy
    return pipeline_cost_dict, path_node_dict


def onsite_to_cluster(spatial_point, price_list, technology, ENVI, standard, facility_ls=30, sewer_ls=30, raster_size=30, pop=2.7,
                      charge=0.08, n=0):
    """
    This function is one of the most import functions of SSM, and is used to merge onsite pattern to community pattern.
    Simulate the merging process from onsite pattern to cluster pattern.

    :param sewer_ls:
    :param pop: 户均人口
    :param raster_size: cell size of raster data set, 30 meter
    :param spatial_point: a nested list, each point represents a household.
                          eg [[ID, X, Y, H, Q, M, 0, [id], [envi, level], tech, c, o&m], ..., [], ...], 0 为管底标高。
    :param price_list: 不同规模设施配套管网系统的户均造价。[p1, p2, p3, p4, p5, p6]
    :param technology: with fields ['id', 'name', 'level', 'scale_min', 'scale_max', 'cost', 'o&m']
    :param ENVI: 水功能区划数据，list[[id, x, y, envi], [...], ...]
    :param standard: 地方/国家农村WWTP出水排放标准
    :param facility_ls: lifespan of sewage treatment facility
    :param charge: daily sewage discharge coefficient, m3 per capita, default value is 0.08 m3/day.
    :param n: the index number of element n in nested list of input list.

    :return:
    scenario_list：nested list, [[id, M, facility, DD1, DD2, S_C, S_M, F_C, F_M, C_T], [], ...].  id: 为编号起始值为1；
                   M为标记（H或C），DD为分散度，计算方法不一；facility为设备数量；S_C, S_M, F_C, F_M, C_T分别为设备（管网）造价、运维费用、总价
    cluster_list: nested list [[ID, X, Y, H, Q, M, 0, [id, id2, ...], [envi, level], tech, c, o&m], ..., [], ...].
                  in which H and M need to update
    """
    def s(e):
        """for sort function"""
        return e[7]

    # data pretreatment
    point_dict = list_to_dict(spatial_point, n)
    household_num = len(spatial_point)
    facility_list = deepcopy(spatial_point)
    cluster_dict_o = get_cluster(facility_list)  # the key of scatter points is "L0", only ID
    if 'L0' in cluster_dict_o:
        cluster_num = len(cluster_dict_o) - 1 + len(cluster_dict_o['L0'])  # number of cluster & scatter
    else:
        cluster_num = len(cluster_dict_o)
    cnt, output_list = 1, []

    # 初始  decentralized pattern
    DD1, DD2 = degree_of_dec(household_num, facility_list, charge, pop)
    facility_list = assign_attribute(facility_list, ENVI, raster_size)  # update envi
    facility_list = get_standard(facility_list, standard, charge)  # update level
    facility_list = system_option(facility_list, technology, charge, facility_ls)  # select technology
    facility_c = sum([math.ceil(i[4] * len(i[7]) * charge) * i[10] for i in facility_list])
    facility_o = sum([math.ceil(i[4] * len(i[7]) * charge) * i[11] for i in facility_list])
    sewer_c, sewer_o = pipe_cost_onsite(facility_list, price_list, sewer_ls)
    cost_t = facility_c + facility_o + sewer_c + sewer_o
    output_list.append([cnt, 'H', len(facility_list), DD1, DD2, sewer_c, sewer_o, facility_c, facility_o, cost_t])

    # iteration merging test
    while len(facility_list) > cluster_num:
        cnt += 1
        cluster_dict = get_cluster(facility_list)
        merged_list = []
        # 合并，用于下一轮迭代计算
        sewer_c, sewer_o, facility_c, facility_o = 0, 0, 0, 0
        # handle cluster point
        for k, v in cluster_dict.items():  # 每轮迭代处理一个聚集区
            if k == 'L0':
                # handle scatter point
                scatter_list = cluster_dict['L0']
                onsite_sc, onsite_so = pipe_cost_onsite(scatter_list, price_list, sewer_ls)
                scatter_list = assign_attribute(scatter_list, ENVI, raster_size)  # update envi
                scatter_list = get_standard(scatter_list, standard, charge)  # update level
                scatter_list = system_option(scatter_list, technology, charge, facility_ls)  # select technology
                onsite_fc = sum([math.ceil(i[4] * len(i[7]) * charge) * i[10] for i in scatter_list])
                onsite_fo = sum([math.ceil(i[4] * len(i[7]) * charge) * i[11] for i in scatter_list])
                merged_list.extend(scatter_list)
                sewer_c += onsite_sc
                sewer_o += onsite_so
                facility_c += onsite_fc
                facility_o += onsite_fo

            elif k != 'L0':
                to_merge_point = v  # household in a residential area
                nearest_point = get_merge_point(to_merge_point)
                facility = to_one_point(nearest_point, point_dict)
                facility = assign_attribute(facility, ENVI, raster_size)  # update envi
                facility = get_standard(facility, standard, charge)  # update level
                facility = system_option(facility, technology, charge, facility_ls)
                f_c = sum([math.ceil(i[4] * len(i[7]) * charge) * i[10] for i in facility])
                f_o = sum([math.ceil(i[4] * len(i[7]) * charge) * i[11] for i in facility])
                s_c, s_o = pipe_cost_onsite(facility, price_list, sewer_ls)
                sewer_c += s_c
                sewer_o += s_o
                facility_c += f_c
                facility_o += f_o
                merged_list.extend(facility)

        facility_num = len(merged_list)
        dd1, dd2 = degree_of_dec(household_num, merged_list, charge, pop)   # check DD1, todo
        total_cost = sewer_c + sewer_o + facility_c + facility_o
        output_list.append([cnt, 'O', facility_num, dd1, dd2, sewer_c, sewer_o, facility_c, facility_o, total_cost])
        facility_list = merged_list

        # 遍历每个聚集区，导致目标区域内cluster_list编号有重复
        facility_list.sort(key=s)
        CNT = 500000
        for c in facility_list:
            if c[5] != 'L0':
                CNT += 1
                ID_M = "M%d" % CNT
                c[0] = ID_M
            else:
                continue
    return output_list, facility_list


def simulate_to_centralized(cluster_point, cluster_dict):
    """
    This function is used to generate the possible scenarios from cluster pattern to centralized pattern.
    通过获取未处理节点的最邻近节点，

    :param cluster_dict: list_to_dict(household_point)
    :param cluster_point: nested point, eg [[ID, X, Y, H, Q, M, 0,[id1, id2, id3], [envi, level], tech, c, o&m], ..., [], ...]

    :return:
    scenario_list: nested list, each element (list) stored ID of the cluster/scatter point which served by one facility
                   [[ID1, ID2, ...], [...], ...]
    """
    # 数据预处理
    scenario_list = []
    cluster_dict0 = deepcopy(cluster_dict)
    scenario = deepcopy(cluster_point)
    # 遍历cluster node计算nearest point
    while len(scenario) > 1:
        nearest_point_list = get_merge_point(scenario)
        scenario = to_one_point(nearest_point_list, cluster_dict0)
        scenario_list.append([s[7] for s in scenario])
    return scenario_list


def to_centralized(point_dict, cluster_point, DEM, type_list, price_list, sewer_cost, technology, ENVI, st, coeff=1.5,
                   depth_min=0.5, wide=0.5, depth_max=4, slope=0.005, lifespan=30, f_ls=30, dis=50, cellsize=30, pop=2.7,
                   charge=0.08, n=0):
    """
    This function is used to test the cost variation for RuST from cluster pattern to centralized pattern by simulating
    the possible pattern of RuST use a merging algorithm. The key input of this function is the cluster point represents
    cluster pattern, and abort algorithm when the centralized pattern is yield.

    :param point_dict: list_to_dict(household_point)
    :param cluster_point: nested list, each point represents a cluster or scatter point.
                          [[ID, X, Y, H, Q, M, T, [id, id2, ...], [E, L], tech, c, o&m], ..., [], ...].
    :param DEM: 由输入DEM数据读取的嵌套列表数据。
    :param type_list: list (ascending order) with standard specifications of pipe, (mm)
    :param price_list: price of earthwork, pipe, drop well, inspection well, pump,eg[12,{type: 59, ...}, 32, 44, 2344]
    :param sewer_cost: sewer cost of cluster pattern. [sewer_c, sewer_o]
    :param technology: with fields ['id', 'name', 'level', 'scale_min', 'scale_max', 'cost', 'o&m']
    :param ENVI: 水功能区划数据，nested list[[id, x, y, envi], [...], ...]
    :param st: 地方/国家农村WWTP出水排放标准
    :param coeff: used to calculate the scanning distance of the given point P, default value is 1.5.
    :param depth_min: the minimum depth of pipeline, the default value is 0.5 meter.
    :param wide: wide of trench, meter
    :param depth_max: the maximum depth of pipeline, the default value is 4 meter.
    :param slope: the slope of pipeline and the default value is 0.005. 管道设计坡度
    :param lifespan: life span of pipeline, default value is 30 year.
    :param f_ls: lifespan of sewage treatment facility
    :param dis: distance between two inspection well, default value is 50 meter.
    :param cellsize: cell size of envi data set, default value is 30 meter.
    :param pop: 户均人口
    :param charge: daily sewage discharge coefficient, m3 per capita, default value is 0.08 m3/day.
    :param n: n: the index number of element n in nested list of input list.

    :return:
    output_list: nested list, [[id, M, facility, DD1, DD2, S_C, S_M, F_C, F_M, C_T], [], ...].  id: 为编号起始值为1；
                 M为标记（H或C），DD为分散度，计算方法不一；facility为设备数量；S_C, S_M, F_C, F_M, C_T分别为设备（管网）造价、
                 运维费用、总价.
    """
    # data pretreatment
    cnt = len(cluster_point)
    # simulating from cluster pattern to clustered pattern
    scenario_list = simulate_to_centralized(cluster_point, point_dict)  # all patterns between cluster pattern and centralized pattern

    cluster_sewer_c, cluster_sewer_o = sewer_cost  # 聚集区内部的管网管网设施造价
    household = 0
    output_list = []
    for x in cluster_point:  # get the number of households
        household += len(x[7])
    # Raise error
    if household == 0:
        raise Exception("ValueError! Number of Household = 0, please check cluster_point list")

    # iteration scenario.
    for scenario in scenario_list:  # 存储元素为ID号，非完整节点信息
        scenario_cluster = []
        facility_num = len(scenario)
        # 计算管网系统造价
        for cluster in scenario:  # 根据节点ID，获取详细信息。cluster 为一个WWTP所服务的对象。
            node_inf_list = [point_dict[ID] for ID in cluster]
            scenario_cluster.append(node_inf_list)  # households served by one
        path_list = get_sewer_system(scenario_cluster, DEM, type_list, price_list, coeff, depth_min, wide, depth_max,
                                     slope, cellsize, lifespan, dis)  # path_dict   path_list

        pipeline_cost_dict, path_node_dict = get_collection_system_cost(path_list, DEM, cluster_point, type_list,
                                                                        price_list, charge, depth_min, wide, depth_max,
                                                                        slope, lifespan)

        pipe_cost_c = sum([c[1] for c in pipeline_cost_dict.values()]) + cluster_sewer_c
        pipe_cost_o = sum([c[2] for c in pipeline_cost_dict.values()]) + cluster_sewer_o
        pipe_cost_t = pipe_cost_c + pipe_cost_o

        # 计算分散度
        merge_point = []
        for c in scenario:   # household id which served by one WWTP
            cluster_households = c[:]
            merge_point.append(['ID', 0, 0, 0, pop, 'M', 0, cluster_households])
        DD1, DD2 = degree_of_dec(household, merge_point, charge, pop)

        # calculate facility cost
        _, source_nodes = get_sewer_path(path_list)  # source_nodes = {WWTP1: [ID1, ID2, ...], ...}
        cluster_dict = list_to_dict(cluster_point, n)
        WWTP_list = []
        for WWTP, source in source_nodes.items():
            wwtp = cluster_dict[WWTP]
            id_7 = []
            for ID in source:
                id_7.extend(cluster_dict[ID][7])
            wwtp[7] = list(set(id_7))
            WWTP_list.append(wwtp)
        WWTP_list = assign_attribute(WWTP_list, ENVI, cellsize)
        WWTP_list = get_standard(WWTP_list, st, charge)
        WWTP_list = system_option(WWTP_list, technology, charge, f_ls)
        facility_c = sum([math.ceil(w[4] * len(w[7]) * charge) * w[10] for w in WWTP_list])
        facility_o = sum([math.ceil(w[4] * len(w[7]) * charge) * w[11] for w in WWTP_list])
        facility_t = facility_c + facility_o
        output_list.append([cnt, 'C', facility_num, DD1, DD2, pipe_cost_c, pipe_cost_o, facility_c, facility_o,
                            pipe_cost_t + facility_t])
        cnt += 1
        gc.collect()
    return output_list


# ----------------------------------------------------------------------------------------------------------------------
# write out
# ----------------------------------------------------------------------------------------------------------------------


def get_output_dict(output_list, household):
    """
    将输出列表转化为字典，以便于输出到文件、或图片输出。

    :param household: number of household which served by sewage treatment facility within the considered area.
    :param output_list: output data of SSM,  [[CNT, 'C', NUM, dd1, DD2, pipe_c, pipe_o, wwtp_c, wwtp_o, total], [], ...]

    :return:
    output_dict: 以CNT, 'C', NUM, dd1, DD2, pipe_c, pipe_o, wwtp_c, wwtp_o, total等为字典的键，便于构造DataFrame
    """
    CNT = [o[0] for o in output_list]
    Marker = [o[1] for o in output_list]
    Number_facility = [o[2] for o in output_list]
    DD_new = [o[3] for o in output_list]
    DD = [o[4] for o in output_list]
    pipe_c = [o[5] / household for o in output_list]
    pipe_o = [o[6] / household for o in output_list]
    wwtp_c = [o[7] / household for o in output_list]
    wwtp_o = [o[8] / household for o in output_list]
    sewer_cost = [(o[5] + o[6]) / household for o in output_list]
    WWTP_cost = [(o[7] + o[8]) / household for o in output_list]
    total = [o[9] / household for o in output_list]
    output_dict = {"01_CNT": CNT, "02_Marker": Marker, "03_Number_facility": Number_facility, "04_DD_new": DD_new,
                   "05_DD": DD, "06_Sewer_con": pipe_c, "07_Sewer_o&m": pipe_o, "08_wwtp_con": wwtp_c,
                   "09_wwtp_o&m": wwtp_o, "10_Sewer_cost": sewer_cost, "11_WWTP_cost": WWTP_cost,
                   "12_Total_cost": total}
    return output_dict


def write_to_file(output_dict, output_file):
    """
    This function is used to write output to a excel file.

    :param output_dict:
    :param output_file: 输出文件及保存目录。

    :return:
    output_file: file with detail information of each simulate scenario.
    """
    output = pd.DataFrame(output_dict)
    output.to_excel(output_file, header=True, sheet_name='SSM_Results')
    return


if __name__ == '__main__':
    a1 = {'D10000347': ['D10000347', 'M10004'], 'D10000346': ['D10000346', 'M10004'], 'D10000348': ['D10000348', 'D10000376'], 
    'M10004': ['M10004', None], 'D10000318': ['D10000318', 'M10004'], 
    'H100034': ['H100034', 'D10000348'], 'D10000376': ['D10000376', 'D10000347']}
    s1 = 'M10004'
    e1 = 'H100034'
    p = get_path(s1, e1, a1)
    print(p)
    a2 = {'D10000891': ['D10000891', 'H100092'], 'H100092': ['H100092', None], 'M10005': ['M10005', 'D10000891']}
    s2 = 'H100092'
    e2 = 'M10005'
    p2 = get_path(s2, e2, a2)
    print(p2)
