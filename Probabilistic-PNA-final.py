import copy
import datetime
from collections import defaultdict
import pandas as pd
import numpy as np


def read_graph_info():
    startpts_path = '.../startpts.txt'
    endpts_path = '.../endpts.txt'
    nodes_path = '.../nodes.txt'
    midpts_path = '.../mid_points.txt'

    f = open(nodes_path, 'r', encoding='UTF-8')
    text = f.readlines()
    nodes = []
    for i in range(1, len(text)):
        nodes.append([round(float(text[i].split(",")[-2]), 6), round(float(text[i].split(",")[-1]), 6)])
    f = open(startpts_path, 'r', encoding='UTF-8')
    text1 = f.readlines()
    f = open(endpts_path, 'r', encoding='UTF-8')
    text2 = f.readlines()
    f = open(midpts_path, 'r', encoding='UTF-8')
    text3 = f.readlines()

    startpts = []
    endpts = []
    midpts = []
    edges_info = []
    for i in range(1, len(text3)):
        startpts.append([round(float(text1[i].split(",")[-2]), 6), round(float(text1[i].split(",")[-1]), 6)])
        endpts.append([round(float(text2[i].split(",")[-2]), 6), round(float(text2[i].split(",")[-1]), 6)])
        midpts.append([round(float(text3[i].split(",")[-2]), 6), round(float(text3[i].split(",")[-1]), 6)])
        # timecost of edge,oneway
        edges_info.append([round(float(text3[i].split(",")[4]), 3), text3[i].split(",")[5]])

    return nodes, startpts, midpts, endpts, edges_info


# 生成有向图的邻接表，需要考虑oneway字段
# Generate the adjacency table of a directed graph, need to consider the oneway field.
def generate_adjacent_list(nodes, startpts, endpts, edges_info):
    # 空间连接里包含的拓扑关系未考虑单向边。
    # The topological relationship contained in the spatial connection does not consider unidirectional edges.
    spatial_join_path = '.../空间连接.txt'
    f = open(spatial_join_path, 'r', encoding='UTF-8')
    text = f.readlines()
    print(len(text))
    print(text[0].split(","))
    # node_fid
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
                if edge_fid[i] == edge and node_fid[i] != edges_list.index(edges):
                    init_graph[edges_list.index(edges)].append(
                        [node_fid[i], edge_fid[i]])

    def dist(pt1, pt2):
        return 1 if np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2) < 1 else 0

    graph = [[] for i in range(len(nodes))]
    for i in range(len(graph)):
        sub_graph = []
        for j in range(len(startpts)):
            if edges_info[j][1] == "1":
                if dist(startpts[j], nodes[i]) or dist(endpts[j], nodes[i]):
                    for k in range(len(init_graph[i])):
                        if j == init_graph[i][k][1]:
                            # 邻接点，邻接边，时间成本(Adjacency vertex, adjacency edge, time cost for adjacency edge)
                            sub_graph.append(
                                [init_graph[i][k][0], init_graph[i][k][1], edges_info[j][0]])
            if edges_info[j][1] == "FT":
                if dist(startpts[j], nodes[i]):
                    for k in range(len(init_graph[i])):
                        if j == init_graph[i][k][1]:
                            # 邻接点，邻接边，时间成本(Adjacency vertex, adjacency edge, time cost for adjacency edge)
                            sub_graph.append(
                                [init_graph[i][k][0], init_graph[i][k][1], edges_info[j][0]])
            if edges_info[j][1] == "TF":
                if dist(endpts[j], nodes[i]):
                    for k in range(len(init_graph[i])):
                        if j == init_graph[i][k][1]:
                            # 邻接点，邻接边，时间成本(Adjacency vertex, adjacency edge, time cost for adjacency edge)
                            sub_graph.append(
                                [init_graph[i][k][0], init_graph[i][k][1], edges_info[j][0]])
        # 输出被去掉的伪邻接节点(Output the removed pseudo adjacency node)
        # sub_list1 = edges_list[i]
        # sub_list2 = list(list(zip(*sub_graph))[1])
        # print(i, list(set(sub_list1) - set(sub_list2)))
        graph[i] = sub_graph
    return graph, len(edges_info)


def cal_angle_cos(v1, v2):
    num = float(np.dot(v1, v2))
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    return num / denom


def hx(pt1, pt2):
    # return (abs(pt1[0] - pt2[0]) + abs(pt1[1] - pt2[1])) / 20
    return np.sqrt((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2) / 20


def Astar_with_turn_cost(graph, nodes_pos, start_node, end_node):
    # 待判断节点集合(节点索引，起点到节点最短路径节点集，起点到节点最短路径边权集，起点到该节点最小时间成本g(x)，
    # 该节点到终点预估时间成本h(x)和g(x)两者之和是f(x)
    # Vertex set to be judged:openlist(vertex index, vertex set of the shortest path from the starting point to the vertex,
    # edge weight set of the shortest path from the starting point to the vertex,
    # minimum time cost g(x) from the starting point to the vertex,
    # The sum of the estimated time cost h(x) from the vertex to the end point and g(x) is f(x)
    openlist = [[start_node, [], [], 0, hx(nodes_pos[start_node], nodes_pos[end_node])]]
    # 已经判断节点集合(Vertex set determined:openlist)
    closelist = []
    # 直到openlist为空或者nmin等于终点时循环结束
    # The loop ends until openlist is empty or nmin[0] is equal to the end point
    while openlist and openlist[0][0] != end_node:
        # 查找openlist中f(x)最小的节点，记为nmin(Find the node with the smallest f(x) in openlist and record it as nmin)
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
        # 查找Nmin的所有不属于Closelist的邻居节点(Find all neighbor nodes of nmin that do not belong to closelist)
        closelist_nodes = list(list(zip(*list(closelist)))[0])
        for i in range(len(graph[nmin[0]])):
            if graph[nmin[0]][i][0] not in closelist_nodes:
                adjacent_nodes.append([graph[nmin[0]][i][0], graph[nmin[0]][i][2]])
        # 如果邻接节点不为空(If adjacent vertices is not empty)
        if adjacent_nodes:
            for i in range(len(adjacent_nodes)):
                # 按照原始A*算法逻辑，加入转弯成本后会出现计算错误的情况，故对原始算法进行改进，
                # 不对比较长的子路径进行更改，因为比较长的子路径可能恰好是最短路径的一部分。
                # According to the logic of the original a * algorithm,
                # there will be calculation errors after adding the turning cost of vertex.
                # Therefore, the original algorithm is improved and the longer sub path is not changed,
                # because the longer sub path may just be part of the shortest path.
                turn_cost = 0
                if nmin[0] != start_node:
                    vec1 = [nodes_pos[nmin[0]][0] - nodes_pos[nmin[1][-1]][0],
                            nodes_pos[nmin[0]][1] - nodes_pos[nmin[1][-1]][1]]
                    vec2 = [nodes_pos[adjacent_nodes[i][0]][0] - nodes_pos[nmin[0]][0],
                            nodes_pos[adjacent_nodes[i][0]][1] - nodes_pos[nmin[0]][1]]
                    turn_cost = cal_turn_cost(vec1, vec2)
                trace_nodes = nmin[1] + [nmin[0]]
                gx = nmin[3] + adjacent_nodes[i][1] + turn_cost
                trace_costs = nmin[2] + [gx]
                fx = gx + hx(nodes_pos[adjacent_nodes[i][0]], nodes_pos[end_node])
                openlist.append([adjacent_nodes[i][0], trace_nodes, trace_costs, gx, fx])
        openlist.sort(key=lambda x: x[4])

    # 输出最短路径点集合和最小时间成本（Output:shortest time,shortest path vertices,time cost list）
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
            # Vertex, another adjacent vertex of adjacent edge, adjacent edge of vertex, minimum time cost of adjacent edge
            all_weight_edges.append((i, graph[i][j][0], graph[i][j][1], graph[i][j][2]))

    # 最小时间成本(minimum time cost:2c)
    c = Astar_with_turn_cost(graph, nodes, start_node, end_node)[0] / 2
    # 最大时间预算(maximum time budget:2a)
    a = multiple * c
    print("c,a", c, a)
    print("optimal_path", Astar_with_turn_cost(graph, nodes, start_node, end_node))

    # 若某节点k的time_sk和time_ke之和大于2a,则无法通过任何边抵达。
    # If the sum of time_sk and time_ke of a vertex k is greater than 2a, it cannot be reached through any edge.
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

    node_time_info = dict(zip([i for i in range(len(nodes))], [[] for _ in range(len(nodes))]))

    for i in range(len(all_weight_edges)):
        if all_weight_edges[i][0] not in node_not_in_ppn and all_weight_edges[i][1] not in node_not_in_ppn:
            modified_graph = modify_graph(all_weight_edges[i][0], all_weight_edges[i][2])
            # 临时更新nodes坐标集(Temporarily update the vertices coordinate set)
            nodes.append(midpts[all_weight_edges[i][2]])

            # 由于转弯成本以边向量构建，边截断后，计算的转弯成本可能会变，以不截断的为准（截断后向量变化，角度变化，右转可能变成了直行）
            # Since the turning cost is constructed by the edge vector, after the edge is truncated,
            # the calculated turning cost may change, and the one that is not truncated shall prevail
            # (the vector changes after the truncation, the angle changes, and the right turn may become straight)
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

                origin_turncost1 = cal_turn_cost(vec1, vec2)
                vec1 = [nodes[node2][0] - nodes[node1][0], nodes[node2][1] - nodes[node1][1]]
                vec2 = [nodes[node4][0] - nodes[node2][0], nodes[node4][1] - nodes[node2][1]]

                new_turncost1 = cal_turn_cost(vec1, vec2)
                # 新减去旧，原来的时间减去这个值（新成本大于旧成本，减去多出来的新城本，反之加上多减去的旧成本）
                # new turning cost - origin turning cost
                time_sk -= new_turncost1 - origin_turncost1

            if len(time_ke_path) > 2:
                node1 = time_sk_path[-2]
                node2 = time_ke_path[1]
                node3 = time_ke_path[2]
                node4 = time_ke_path[0]
                # 2-1x3-2 2-4x3-2
                vec1 = [nodes[node2][0] - nodes[node1][0], nodes[node2][1] - nodes[node1][1]]
                vec2 = [nodes[node3][0] - nodes[node2][0], nodes[node3][1] - nodes[node2][1]]

                origin_turncost2 = cal_turn_cost(vec1, vec2)
                vec1 = [nodes[node2][0] - nodes[node4][0], nodes[node2][1] - nodes[node4][1]]
                vec2 = [nodes[node3][0] - nodes[node2][0], nodes[node3][1] - nodes[node2][1]]

                new_turncost2 = cal_turn_cost(vec1, vec2)
                # 新减去旧，原来的时间减去这个值（新成本大于旧成本，减去多出来的新城本，反之加上多减去的旧成本）
                # new turning cost - origin turning cost
                time_ke -= new_turncost2 - origin_turncost2

            # 由于k是路段中点无转弯成本，则有tsk + tke = tse
            # Since k is the midpoint of the road section without turning cost, there is tsk+tke = tse.
            total_time = time_sk + time_ke
            # 恢复nodes坐标集(Restore nodes coordinate set)
            nodes.pop(-1)

            # 由total_time(tse)计算s到进入路段的节点k的时间成本modified_time_sk（即下面的tsk）和出路段的节点kmodified_time_ke（即下面的tke）
            # By total_time(tse), calculates the time cost from s to the vertex k(entry road section):modified_time_sk (i.e. tsk below)
            # andkthe time cost from k(exit road section) to e: modified_time_sk (i.e. tke below)
            # 这里和直接计算tsk和tke的结果存在细微差别(There is a slight difference here and the result of directly computing tsk and tke)
            if total_time <= 2 * a:
                if all_weight_edges[i][0] == time_sk_path[-2]:
                    modified_time_sk = time_sk - edges_info[all_weight_edges[i][2]][0] / 2
                    modified_time_ke = time_ke + edges_info[all_weight_edges[i][2]][0] / 2
                    if len(time_sk_path) > 2:
                        vec1 = [nodes[time_sk_path[-2]][0] - nodes[time_sk_path[-3]][0],
                                nodes[time_sk_path[-2]][1] - nodes[time_sk_path[-3]][1]]
                        vec2 = [nodes[time_ke_path[1]][0] - nodes[time_sk_path[-2]][0],
                                nodes[time_ke_path[1]][1] - nodes[time_sk_path[-2]][1]]
                        turn_cost = cal_turn_cost(vec1, vec2)
                        modified_time_sk -= turn_cost
                    if not node_time_info[all_weight_edges[i][0]]:
                        node_time_info[all_weight_edges[i][0]] = [modified_time_sk, modified_time_ke, total_time]
                    else:
                        if total_time <= node_time_info[all_weight_edges[i][0]][2]:
                            node_time_info[all_weight_edges[i][0]] = [modified_time_sk, modified_time_ke, total_time]
                if all_weight_edges[i][1] == time_ke_path[1]:
                    modified_time_sk = time_sk + edges_info[all_weight_edges[i][2]][0] / 2
                    modified_time_ke = time_ke - edges_info[all_weight_edges[i][2]][0] / 2
                    if len(time_ke_path) > 2:
                        vec1 = [nodes[time_ke_path[1]][0] - nodes[time_sk_path[-2]][0],
                                nodes[time_ke_path[1]][1] - nodes[time_sk_path[-2]][1]]
                        vec2 = [nodes[time_ke_path[2]][0] - nodes[time_ke_path[1]][0],
                                nodes[time_ke_path[2]][1] - nodes[time_ke_path[1]][1]]
                        turn_cost = cal_turn_cost(vec1, vec2)
                        modified_time_ke -= turn_cost
                    if not node_time_info[all_weight_edges[i][1]]:
                        node_time_info[all_weight_edges[i][1]] = [modified_time_sk, modified_time_ke, total_time]
                    else:
                        if total_time <= node_time_info[all_weight_edges[i][1]][2]:
                            node_time_info[all_weight_edges[i][1]] = [modified_time_sk, modified_time_ke, total_time]

    # 删去字典值为空项(Delete empty dictionary value)
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
        # 这里不一定满足tsk + tke = tse（tsk + tke = tse is not necessarily satisfied here）
        # tsk和tke在tse基础上派生而来（tsk and tke are derived from tse）
        # equation 3
        wse = np.power(np.e, -beta1 * ((tse - 2 * c) / (2 * a - 2 * c)))
        # equation 4
        weight1 = np.power(np.e, -beta2 * (tsk / (a + c)))
        weight2 = np.power(np.e, -beta2 * (tke / (a + c)))
        wk = weight1 + weight2

        # TGDE（i.e. ridge model）
        # wk = 1

        weight = wse * wk
        node_weight_info.append(
            [e, tsk, tke, tse, wse, wk, weight])

    return c, a, node_weight_info


def modify_graph(nodei_index, edge_index):
    # m无转弯成本，从s出发经过m再到e的成本为tsm+tme
    # M has no turning cost. The cost of starting from s, passing through m and then arriving at e is tsm + tme
    new_graph = copy.deepcopy(graph)
    # 将中点m所在边截断成2段（Cut the edge of the midpoint m into 2 sections(eim & emj)）
    # 以i到m到j为例，图新加一项，存放m与j的邻接关系(j，emj，tmj（相比以前减半）)
    # Taking i to m to j as an example, a new element is added to the graph,
    # to store the adjacency relationship between m and j (j,emj,tmj(halved compared with the previous))）
    # 图中索引为i的节点，去掉i与j的邻接关系，增加i与m的邻接关系(m,eim,tim)
    # For the vertex with index i in the graph,
    # remove the adjacency relationship between i and j and increase the adjacency relationship between i and m (m,eim,tim)
    # 新增2个边(add 2 egdes)
    new_edge_index = [len(midpts), len(midpts) + 1]
    # mid_pt:m
    new_node_index = len(nodes)
    # 中点m的新子图，只有1条到j的边(The new subgraph of midpoint m has only 1 edge to j)
    new_sub_graph = []
    for i in range(len(new_graph)):
        if i == nodei_index:
            adjacent_edges = list(list(zip(*new_graph[i]))[1])
            # j
            nodej_index = new_graph[i][adjacent_edges.index(edge_index)][0]
            timecost = new_graph[i][adjacent_edges.index(edge_index)][2]
            # 去掉i与j的邻接关系(Remove the adjacency relationship between i and j)
            new_graph[i].pop(adjacent_edges.index(edge_index))
            # 增加i与m的邻接关系(Increase the adjacency relationship between i and m)
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
    # 权重转概率（Weight to probability）
    # equation 5
    new_weights = []
    for i in range(len(node_weight_info)):
        new_weights.append(node_weight_info[i][-1] / all_sum_weights)

    df = pd.DataFrame({'node': list(list(zip(*node_weight_info))[0]),
                       'tsk': list(list(zip(*node_weight_info))[1]),
                       'tke': list(list(zip(*node_weight_info))[2]),
                       'tse': list(list(zip(*node_weight_info))[3]),
                       'wse': list(list(zip(*node_weight_info))[4]),
                       'wk': list(list(zip(*node_weight_info))[5]),
                       'weight': list(list(zip(*node_weight_info))[6]),
                       'probability': new_weights})
    save_path = "c:/Users/RTX2080Ti/Desktop/sheets/点权表" + str(start_node) + "-" + str(end_node) + ".xlsx"
    df.to_excel(save_path)


def cal_turn_cost(vec1, vec2):
    vec_product = vec2[0] * vec1[1] - vec1[0] * vec2[1]
    num = float(np.dot(vec2, vec1))  # 向量点乘
    denom = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    angle = 180 / np.pi * np.arccos(1 if abs(num / denom) >= 1 else num / denom)
    # 直行或者掉头（Go straight or U-turn）
    if angle < 45 or angle > 135:
        turn_cost = 0 if num > 0 else 40
    # 左转或者右转（Turn left or right）
    else:
        turn_cost = 3 if vec_product > 0 else 15
    return turn_cost


if __name__ == '__main__':
    start_t = datetime.datetime.now()
    nodes, startpts, midpts, endpts, edges_info = read_graph_info()
    graph, edge_count = generate_adjacent_list(nodes, startpts, endpts, edges_info)

    # anchor point id
    start_node = 86
    end_node = 81
    # a = 1.5c
    multiple = 1.5

    c, a, node_weight_info = generate_weight()
    modify_node_weights()
    end_t = datetime.datetime.now()
    sec = (end_t - start_t).total_seconds()
    print("所用时间", sec)
    exit()
