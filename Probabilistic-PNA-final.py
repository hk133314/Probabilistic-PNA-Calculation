import copy
import datetime
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import numpy as np


def read_txt():
    startpts_path = 'c:/Users/RTX2080Ti/Desktop/西安二环内路网1/startpts1.txt'
    endpts_path = 'c:/Users/RTX2080Ti/Desktop/西安二环内路网1/endpts1.txt'
    nodes_path = 'c:/Users/RTX2080Ti/Desktop/西安二环内路网1/nodes1.txt'
    midpts_path = 'c:/Users/RTX2080Ti/Desktop/西安二环内路网1/mid_points2.txt'

    f = open(nodes_path, 'r', encoding='UTF-8')
    text = f.readlines()
    print(len(text))
    print(text[0].split(","))

    nodes = []
    for i in range(1, len(text)):
        nodes.append([round(float(text[i].split(",")[-2]), 6), round(float(text[i].split(",")[-1]), 6)])

    f = open(startpts_path, 'r', encoding='UTF-8')
    text1 = f.readlines()
    print(len(text1))
    print(text1[0].split(","))

    f = open(endpts_path, 'r', encoding='UTF-8')
    text2 = f.readlines()
    print(len(text2))
    print(text2[0].split(","))

    f = open(midpts_path, 'r', encoding='UTF-8')
    text3 = f.readlines()
    print(len(text3))
    print(text3[0].split(","))

    startpts = []
    endpts = []
    midpts = []
    edges_info = []
    for i in range(1, len(text3)):
        startpts.append([round(float(text1[i].split(",")[-2]), 6), round(float(text1[i].split(",")[-1]), 6)])
        endpts.append([round(float(text2[i].split(",")[-2]), 6), round(float(text2[i].split(",")[-1]), 6)])
        midpts.append([round(float(text3[i].split(",")[-2]), 6), round(float(text3[i].split(",")[-1]), 6)])
        # timecost,oneway
        edges_info.append([round(float(text3[i].split(",")[4]), 3), text3[i].split(",")[5]])

    return nodes, startpts, midpts, endpts, edges_info


# 有向图的邻接表，需要考虑oneway字段。
def generate_adjacent_list(nodes, startpts, endpts, edges_info):
    # 空间连接里包含的拓扑关系未考虑单向边。
    spatial_join_path = 'c:/Users/RTX2080Ti/Desktop/西安二环内路网1/空间连接1.txt'
    f = open(spatial_join_path, 'r', encoding='UTF-8')
    text = f.readlines()
    print(len(text))
    print(text[0].split(","))
    # node_fid是节点的索引
    node_fid = []
    edge_fid = []
    # node_fid_list:target_fid edge_fid_list:join_fid edge_info:weight
    for i in range(1, len(text)):
        node_fid.append(int(text[i].split(",")[2]))
        edge_fid.append(int(text[i].split(",")[3]))

    nodes_list = defaultdict(list)
    for i, e in enumerate(node_fid):
        nodes_list[e].append(e)
    nodes_list = list(nodes_list.values())

    edges_list = [[] for i in range(len(nodes_list))]
    count = 0
    for i in range(len(nodes_list)):
        for j in range(len(nodes_list[i])):
            edges_list[i].append(edge_fid[count])
            count += 1

    init_graph = [[] for i in range(len(nodes_list))]
    for edges in edges_list:
        for edge in edges:
            for i in range(len(edge_fid)):
                # 注意node_fid是从1开始的还是从0开始的
                if edge_fid[i] == edge and node_fid[i] != edges_list.index(edges):
                    init_graph[edges_list.index(edges)].append(
                        [node_fid[i], edge_fid[i]])

    def dist(pt1, pt2):
        return 1 if np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2) < 1 else 0

    graph = [[] for i in range(len(nodes))]
    for i in range(len(graph)):
        sub_graph = []
        # print(init_graph[i])
        for j in range(len(startpts)):
            if edges_info[j][1] == "1":
                if dist(startpts[j], nodes[i]) or dist(endpts[j], nodes[i]):
                    # print(i, j)
                    for k in range(len(init_graph[i])):
                        if j == init_graph[i][k][1]:
                            # 邻接点，邻接边，时间成本
                            sub_graph.append(
                                [init_graph[i][k][0], init_graph[i][k][1], edges_info[j][0]])
            if edges_info[j][1] == "FT":
                if dist(startpts[j], nodes[i]):
                    # print(i, j)
                    for k in range(len(init_graph[i])):
                        if j == init_graph[i][k][1]:
                            # 邻接点，邻接边，时间成本
                            sub_graph.append(
                                [init_graph[i][k][0], init_graph[i][k][1], edges_info[j][0]])
            if edges_info[j][1] == "TF":
                if dist(endpts[j], nodes[i]):
                    # print(i, j)
                    for k in range(len(init_graph[i])):
                        if j == init_graph[i][k][1]:
                            # 邻接点，邻接边，时间成本
                            sub_graph.append(
                                [init_graph[i][k][0], init_graph[i][k][1], edges_info[j][0]])
        # 输出被去掉的伪邻接节点
        # sub_list1 = edges_list[i]
        # sub_list2 = list(list(zip(*sub_graph))[1])
        # print(i, list(set(sub_list1) - set(sub_list2)))
        graph[i] = sub_graph
    return graph, len(edges_info)


def cal_angle_cos(v1, v2):
    num = float(np.dot(v1, v2))  # 向量点乘
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    return num / denom


# Astar h(x)：曼哈顿距离（速度选个比较大的，例如20m/s）
def hx(pt1, pt2):
    # return (abs(pt1[0] - pt2[0]) + abs(pt1[1] - pt2[1])) / 30
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2) / 20


def Astar_with_turn_cost(graph, nodes_pos, start_node, end_node):
    # 待判断节点集合(节点索引，起点到节点最短路径节点集，起点到节点最短路径边权集，起点到该节点最小时间成本g(x)，
    # 该节点到终点预估时间成本h(x)和起点到该节点最小时间成本g(x)两者之和f(x)
    openlist = [[start_node, [], [], 0, hx(nodes_pos[start_node], nodes_pos[end_node])]]
    # 已经判断节点集合
    closelist = []
    # 直到openlist为空或者nmin等于终点时循环结束
    while openlist and openlist[0][0] != end_node:
        # 查找openlist中f(x)最小的节点，记为nmin
        nmin = openlist[0]
        if not closelist or nmin[0] not in list(list(zip(*closelist))[0]):
            closelist.append(nmin)
        if closelist and nmin[0] in list(list(zip(*closelist))[0]):
            old_idx = list(list(zip(*closelist))[0]).index(nmin[0])
            old_gx = closelist[old_idx][3]
            if nmin[3] <= old_gx:
                closelist[old_idx] = nmin
        openlist.pop(0)
        adjacent_nodes = []
        # 查找Nmin的所有不属于Closelist的邻居节点
        closelist_nodes = list(list(zip(*list(closelist)))[0])
        for i in range(len(graph[nmin[0]])):
            if graph[nmin[0]][i][0] not in closelist_nodes:
                adjacent_nodes.append([graph[nmin[0]][i][0], graph[nmin[0]][i][2]])
        # 如果邻接节点不为空
        if adjacent_nodes:
            # openlist节点集合
            for i in range(len(adjacent_nodes)):
                # 原始A*算法加入转弯成本后会出现计算错误的情况，故对原始算法进行改进，
                # 不对比较长的子路径进行更改，因为比较长的子路径可能恰好是最短路径的一部分。
                turn_cost = 0
                if nmin[0] != start_node:
                    vec1 = [nodes_pos[nmin[0]][0] - nodes_pos[nmin[1][-1]][0],
                            nodes_pos[nmin[0]][1] - nodes_pos[nmin[1][-1]][1]]
                    vec2 = [nodes_pos[adjacent_nodes[i][0]][0] - nodes_pos[nmin[0]][0],
                            nodes_pos[adjacent_nodes[i][0]][1] - nodes_pos[nmin[0]][1]]
                    turn_cost = turn_penalty(vec1, vec2)

                trace_nodes = nmin[1] + [nmin[0]]
                gx = nmin[3] + adjacent_nodes[i][1] + turn_cost
                trace_costs = nmin[2] + [gx]
                fx = gx + hx(nodes_pos[adjacent_nodes[i][0]], nodes_pos[end_node])
                openlist.append([adjacent_nodes[i][0], trace_nodes, trace_costs, gx, fx])
        # 按照gx排序比按照fx排序慢不少，但是效果是一样的
        openlist.sort(key=lambda x: x[4])

    # 输出最短路径点集合和最小时间成本
    if openlist and openlist[0][0] == end_node:
        closelist.sort(key=lambda x: x[0])
        optimal_time = openlist[0][3]
        optimal_path_nodes = openlist[0][1] + [end_node]
        time_costs = [0] + openlist[0][2]
        old_time_costs = copy.deepcopy(time_costs)
        if len(time_costs) > 1:
            for i in range(1, len(time_costs)):
                time_costs[i] = round(old_time_costs[i] - old_time_costs[i - 1], 3)
        return optimal_time, optimal_path_nodes, time_costs
    if not openlist:
        return -1, [], []


def generate_weight():
    all_weight_edges = []
    for i in range(len(graph)):
        for j in range(len(graph[i])):
            # 节点，边的邻接节点，邻接边，最小时间成本
            all_weight_edges.append((i, graph[i][j][0], graph[i][j][1], graph[i][j][2]))

    # 最小时间成本
    c = Astar_with_turn_cost(graph, nodes, start_node, end_node)[0] / 2
    a = multiple * c
    print("c,a", c, a)
    print("optimal_path", Astar_with_turn_cost(graph, nodes, start_node, end_node))

    # 若某节点的time_sk和time_ke之和大于2a,则无法通过任何边抵达。
    node_not_in_ppn = []
    for i in range(len(nodes)):
        # print(i)
        node_time_sk = Astar_with_turn_cost(graph, nodes, start_node, i)[0]
        if node_time_sk > 2 * a:
            node_not_in_ppn.append(i)
        else:
            node_time_ke = Astar_with_turn_cost(graph, nodes, i, end_node)[0]
            if node_time_sk + node_time_ke > 2 * a:
                node_not_in_ppn.append(i)
    # print(len(node_not_in_ppn), node_not_in_ppn)

    node_time_info = dict(zip([i for i in range(len(nodes))], [[] for _ in range(len(nodes))]))
    ppn_edge_index = []

    for i in range(len(all_weight_edges)):
        if all_weight_edges[i][0] not in node_not_in_ppn and all_weight_edges[i][1] not in node_not_in_ppn:
            modified_graph = modify_graph(all_weight_edges[i][0], all_weight_edges[i][2])
            # 临时更新nodes坐标集
            nodes.append(midpts[all_weight_edges[i][2]])

            # 边截断后，计算的时间成本可能会变（影响转弯成本值，截断后根据角度，直行变成了右转
            time_sk, time_sk_path, timecost1 = Astar_with_turn_cost(modified_graph, nodes, start_node,
                                                                    len(nodes) - 1)
            time_ke, time_ke_path, timecost2 = Astar_with_turn_cost(modified_graph, nodes, len(nodes) - 1,
                                                                    end_node)

            if len(time_sk_path) > 2:
                node1 = time_sk_path[-3]
                node2 = time_sk_path[-2]
                node3 = time_ke_path[1]
                node4 = time_ke_path[0]
                # 2-1x3-2 2-1x4-2
                vec1 = [nodes[node2][0] - nodes[node1][0], nodes[node2][1] - nodes[node1][1]]
                vec2 = [nodes[node3][0] - nodes[node2][0], nodes[node3][1] - nodes[node2][1]]

                origin_turncost1 = turn_penalty(vec1, vec2)
                vec1 = [nodes[node2][0] - nodes[node1][0], nodes[node2][1] - nodes[node1][1]]
                vec2 = [nodes[node4][0] - nodes[node2][0], nodes[node4][1] - nodes[node2][1]]

                new_turncost1 = turn_penalty(vec1, vec2)
                # 新减去旧，原来的时间减去这个值（新成本大于旧成本，减去多出来的新城本，反之加上多减去的旧成本）
                time_sk -= new_turncost1 - origin_turncost1

            if len(time_ke_path) > 2:
                node1 = time_sk_path[-2]
                node2 = time_ke_path[1]
                node3 = time_ke_path[2]
                node4 = time_ke_path[0]
                # 2-1x3-2 2-4x3-2
                vec1 = [nodes[node2][0] - nodes[node1][0], nodes[node2][1] - nodes[node1][1]]
                vec2 = [nodes[node3][0] - nodes[node2][0], nodes[node3][1] - nodes[node2][1]]

                origin_turncost2 = turn_penalty(vec1, vec2)
                vec1 = [nodes[node2][0] - nodes[node4][0], nodes[node2][1] - nodes[node4][1]]
                vec2 = [nodes[node3][0] - nodes[node2][0], nodes[node3][1] - nodes[node2][1]]

                new_turncost2 = turn_penalty(vec1, vec2)
                # 新减去旧，原来的时间减去这个值（新成本大于旧成本，减去多出来的新城本，反之加上多减去的旧成本）
                time_ke -= new_turncost2 - origin_turncost2

            # 由于k是路段中点，则有tsk + tke = tse
            total_time = time_sk + time_ke
            # 恢复nodes坐标集
            nodes.pop(-1)

            # 满足该条件的加入，前提：要确保这个total_time是正确的
            if total_time <= 2 * a:
                # time_sk_path或者time_ke_path长度等于2时(按路段中心点计算时至少为2)
                # 目标节点(不论是all_weight_edges[i][0]还是all_weight_edges[i][1])不会位于其中（起点终点除外）
                # 若位于其中，则会位于time_sk_path的倒数第2项或time_ke_path的第2项。
                # 这里好像不对
                if all_weight_edges[i][0] == time_sk_path[-2]:
                    modified_time_sk = time_sk - edges_info[all_weight_edges[i][2]][0] / 2
                    modified_time_ke = time_ke + edges_info[all_weight_edges[i][2]][0] / 2
                    # modified_time_sk需要减去通过路段时间成本的一半，再减去从目标节点过渡到time_ke_path的转弯成本。
                    # time_sk_path[-2]：目标节点
                    # len(time_sk_path) == 2时为起点，无转弯成本
                    if len(time_sk_path) > 2:
                        vec1 = [nodes[time_sk_path[-2]][0] - nodes[time_sk_path[-3]][0],
                                nodes[time_sk_path[-2]][1] - nodes[time_sk_path[-3]][1]]
                        vec2 = [nodes[time_ke_path[1]][0] - nodes[time_sk_path[-2]][0],
                                nodes[time_ke_path[1]][1] - nodes[time_sk_path[-2]][1]]
                        turn_cost = turn_penalty(vec1, vec2)
                        modified_time_sk -= turn_cost
                    if not node_time_info[all_weight_edges[i][0]]:
                        node_time_info[all_weight_edges[i][0]] = [modified_time_sk, modified_time_ke, total_time]
                    else:
                        if total_time <= node_time_info[all_weight_edges[i][0]][2]:
                            node_time_info[all_weight_edges[i][0]] = [modified_time_sk, modified_time_ke, total_time]
                if all_weight_edges[i][1] == time_ke_path[1]:
                    modified_time_sk = time_sk + edges_info[all_weight_edges[i][2]][0] / 2
                    modified_time_ke = time_ke - edges_info[all_weight_edges[i][2]][0] / 2
                    # 需要减去从目标节点过渡到time_ke_path的转弯成本。
                    # time_ke_path[1]：目标节点
                    # len(time_ke_path) == 2时为终点，无转弯成本
                    if len(time_ke_path) > 2:
                        vec1 = [nodes[time_ke_path[1]][0] - nodes[time_sk_path[-2]][0],
                                nodes[time_ke_path[1]][1] - nodes[time_sk_path[-2]][1]]
                        vec2 = [nodes[time_ke_path[2]][0] - nodes[time_ke_path[1]][0],
                                nodes[time_ke_path[2]][1] - nodes[time_ke_path[1]][1]]
                        turn_cost = turn_penalty(vec1, vec2)
                        modified_time_ke -= turn_cost
                    if not node_time_info[all_weight_edges[i][1]]:
                        node_time_info[all_weight_edges[i][1]] = [modified_time_sk, modified_time_ke, total_time]
                    else:
                        if total_time <= node_time_info[all_weight_edges[i][1]][2]:
                            node_time_info[all_weight_edges[i][1]] = [modified_time_sk, modified_time_ke, total_time]

    # 删去字典值为空项
    for i in range(len(node_time_info)):
        if not node_time_info[i]:
            del node_time_info[i]

    node_weight_info = []
    beta1 = 5
    beta2 = 3
    for i, e in enumerate(node_time_info):
        tsk = node_time_info[e][0]
        tke = node_time_info[e][1]
        tse = node_time_info[e][2]
        # 这里可不一定是tsk + tke = tse
        wse = np.power(np.e, -beta1 * ((tse - 2 * c) / (2 * a - 2 * c)))
        weight1 = np.power(np.e, -beta2 * (tsk / (a + c)))
        weight2 = np.power(np.e, -beta2 * (tke / (a + c)))
        wk = weight1 + weight2

        # TGDE
        wk = 1

        # TGDE不考虑wk
        weight = wse * wk
        node_weight_info.append(
            [e, tsk, tke, tse, wse, wk, weight])

    return c, a, node_weight_info


def modify_graph2(graph, del_edge_info):
    # 深拷贝新列表
    new_graph = copy.deepcopy(graph)
    # 对del_edge_info中邻接节点代表的边逐一进行删除
    for i in range(len(del_edge_info)):
        adjacent_nodes = list(list(zip(*new_graph[del_edge_info[i][0]]))[0])
        new_graph[del_edge_info[i][0]].pop(adjacent_nodes.index(del_edge_info[i][1]))
    return new_graph


def modify_graph(nodei_index, edge_index):
    # 深拷贝新列表
    new_graph = copy.deepcopy(graph)
    # 将中点所在边截断成2段。
    # eij变成eim emj
    # 整个图经修改后，加入了2条新边
    # 以i到m到j为例，图新加一项，存放m与j的邻接关系(j，emj，tmj（相比以前减半）
    # 图中索引为i的节点，去掉i与j的邻接关系，增加i与m的邻接关系(m,eim,tim,tvim)
    # mid_pt = midpts[edge_index]
    # 新增2个边
    new_edge_index = [len(midpts), len(midpts) + 1]
    # 中点m
    new_node_index = len(nodes)
    # 中点m的新子图，只有1条到j的边
    new_sub_graph = []
    for i in range(len(new_graph)):
        if i == nodei_index:
            adjacent_edges = list(list(zip(*new_graph[i]))[1])
            # j
            nodej_index = new_graph[i][adjacent_edges.index(edge_index)][0]
            timecost = new_graph[i][adjacent_edges.index(edge_index)][2]
            # 去掉i与j的邻接关系
            new_graph[i].pop(adjacent_edges.index(edge_index))
            # 增加i与m的邻接关系
            new_graph[i].append([new_node_index, new_edge_index[0], timecost / 2])
            new_sub_graph = [[nodej_index, new_edge_index[1], timecost / 2]]
            break
    new_graph.append(new_sub_graph)
    return new_graph


def modify_node_weights():
    # [i, tsk, tke, tse, wse, wk, weight]
    all_sum_weights = 0
    for i in range(len(node_weight_info)):
        all_sum_weights += node_weight_info[i][-1]
    # 权重转概率
    new_weights = []
    for i in range(len(node_weight_info)):
        new_weights.append(node_weight_info[i][-1] / all_sum_weights)

    df = pd.DataFrame({'node': list(list(zip(*node_weight_info))[0]),
                       # 'x': list(list(zip(*node_time_info))[1]),
                       # 'y': list(list(zip(*node_time_info))[2]),
                       'tsk': list(list(zip(*node_weight_info))[1]),
                       'tke': list(list(zip(*node_weight_info))[2]),
                       'tse': list(list(zip(*node_weight_info))[3]),
                       'wse': list(list(zip(*node_weight_info))[4]),
                       'wk': list(list(zip(*node_weight_info))[5]),
                       'weight': list(list(zip(*node_weight_info))[6]),
                       'probability': new_weights})
    save_path = "c:/Users/RTX2080Ti/Desktop/sheets/点权表" + str(start_node) + "-" + str(end_node) + ".xlsx"
    df.to_excel(save_path)


def turn_penalty(vec1, vec2):
    # vec1：插入的上一条边
    # vec2：待插入边
    # 向量积，大于0说明在vec1右边，顺时针夹角小于180度，反之在vec1左边，顺时针夹角大于180度。
    vec_product = vec2[0] * vec1[1] - vec1[0] * vec2[1]
    num = float(np.dot(vec2, vec1))  # 向量点乘
    denom = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    angle = 180 / np.pi * np.arccos(1 if abs(num / denom) >= 1 else num / denom)
    # 直行或者掉头
    turn_cost = 0
    if angle < 45 or angle > 135:
        # print(vec1, vec2, vec_product, angle)
        turn_cost = 0 if num > 0 else 40
    else:
        turn_cost = 3 if vec_product > 0 else 15
    return turn_cost


if __name__ == '__main__':
    start_t = datetime.datetime.now()
    # 思路：
    # 1.根据Arcmap空间连接工具生成的邻接信息表，生成节点的邻接列表（列表的每项按照节点索引升序排列）。
    nodes, startpts, midpts, endpts, edges_info = read_txt()
    graph, edge_count = generate_adjacent_list(nodes, startpts, endpts, edges_info)

    start_node = 86
    end_node = 81
    multiple = 1.5

    c, a, node_weight_info = generate_weight()
    modify_node_weights()
    end_t = datetime.datetime.now()
    sec = (end_t - start_t).total_seconds()
    print("所用时间", sec)
    exit()
