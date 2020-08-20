from numpy import unique
import scipy as sp

def py_neighbours(_np, _ien):
    result_ele = [None] * _np
    result_node = [None] * _np
    for i in range(_np):
        result_ele[i] = []
        result_node[i] = []
        for e in range(len(_ien)):
            vertex = sp.array([_ien[e, 0], _ien[e, 1], _ien[e, 2]])
            mask = sp.ones(len(vertex), dtype=bool)
            for j in range(3):
                if vertex[j] == i:
                    mask[j] = False
                    result_ele[i].append(e)
                    temp_result = vertex[mask,...]
                    flag0 = 0
                    flag1 = 0
                    for k in result_node[i]:
                        if k == temp_result[0]:
                            flag0 = 1
                        if k == temp_result[1]:
                            flag1 = 1
                    if flag0 == 0:
                        result_node[i].append(temp_result[0])
                    if flag1 == 0:
                        result_node[i].append(temp_result[1])
                    break

    return result_ele, result_node

def py_neighbours2(_np, _ien):
    result_ele = [None] * _np
    result_node = [None] * _np
    for i in range(_np):
        result_ele[i] = []
        result_node[i] = []
        for e in range(len(_ien)):
            vertex = [_ien[e, 0], _ien[e, 1], _ien[e, 2]]
            for j in range(3):
                if vertex[j] == i:
                    result_ele[i].append(e)
                    del vertex[j]
                    result_node[i].extend(vertex)
                    break
        result_node[i] = list(unique(result_node[i]))


    return result_ele, result_node
