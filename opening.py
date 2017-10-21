#!/usr/bin/env python3

# T : total number of leaves in the tree
# leaves : indexes of the provided leaves
# return : the indexes of the elements to provide so that the full Merkle Tree can be computed
#   root is (0, 0), leaves are from (height - 1, 0) to (height - 1, T - 1)
def opening(T, leaves):
    res = []
    from math import log2
    height = int(log2(T)) + 1
    currentHeight = height - 2
    numberOfNodes = T
    nodes = [False] * numberOfNodes
    for leaf in leaves:
        nodes[leaf] = True
    currentNode = 0
    while currentHeight >= 0:
        numberOfNodes /= 2
        # print('height : %d, nodes : %d' % (currentHeight, numberOfNodes))
        while currentNode < numberOfNodes:
            leftNode = 2 * currentNode
            rightNode = leftNode + 1
            # print('left : %d, right : %d' % (leftNode, rightNode))
            if nodes[leftNode] and nodes[rightNode]:
                nodes[currentNode] = True
            elif nodes[leftNode] and not nodes[rightNode]:
                res.append((currentHeight + 1, rightNode))
                nodes[currentNode] = True
            elif not nodes[leftNode] and nodes[rightNode]:
                res.append((currentHeight + 1, leftNode))
                nodes[currentNode] = True
            else:
                nodes[currentNode] = False
            currentNode += 1
        currentHeight -= 1
        currentNode = 0
    return res

print(opening(8, [0, 2, 4, 6]))
print(opening(8, [0, 7]))
print(opening(8, [0]))

# with the same arguments, returns the list of the indexes to provide with the leaves if the Merkle Tree is stored as in the 'merkle_tree' function
def openingForOneArray(T, leaves):
    nodes = opening(T, leaves)
    return list(map(lambda node : 2 ** node[0] - 1 + node[1], nodes))

print(openingForOneArray(8, [0, 2, 4, 6]))
print(openingForOneArray(8, [0, 7]))
print(openingForOneArray(8, [0]))

def opening_2(T, leaves):
    res = []
    from math import log2
    height = int(log2(T)) + 1
    currentHeight = height - 2
    numberOfNodes = T
    nodes = sorted(leaves)
    while currentHeight >= 0:
        numberOfNodes /= 2
        nodeIndex = 0
        parentNodes = []
        while nodeIndex < len(nodes):
            parentNodes.append(int(nodes[nodeIndex] / 2))
            if nodes[nodeIndex] % 2 == 0:
                if nodeIndex + 1 < len(nodes) and nodes[nodeIndex + 1] == nodes[nodeIndex] + 1:
                    nodeIndex += 2
                else:
                    res.append((currentHeight + 1, nodes[nodeIndex] + 1))
                    nodeIndex += 1
            else:
                res.append((currentHeight + 1, nodes[nodeIndex] - 1))
                nodeIndex += 1
        currentHeight -= 1
        nodes = parentNodes.copy()
    return res

print(opening_2(8, [0, 2, 4, 6]) == opening(8, [0, 2, 4, 6]))
print(opening_2(8, [0, 7]) == opening(8, [0, 7]))
print(opening_2(8, [0]) == opening(8, [0]))
