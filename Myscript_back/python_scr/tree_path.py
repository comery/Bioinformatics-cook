import copy


class TreeNode:
    def __init__(self, val):
        self.val = val
        self.left, self.right = None, None


def dfs(node, result, tmp=list()):
    if node is None:
        return

    tmp.append(node)
    # 这里需要用拷贝而不是用 = 赋值，也可以遍历赋值
    tmp1 = copy.deepcopy(tmp)

    if node.left is None and node.right is None:
        result.append([i.val for i in tmp])
        return

    if node.left is not None:
        dfs(node.left, result, tmp)
    # 遍历右子树需要带上不同的变量，否则左子树的tmp和右子树的tmp都指向一块内存
    if node.right is not None:
        dfs(node.right, result, tmp1)

def is_odd(dd):
    if dd % 2 == 0:
        return True
    else:
        return False

if __name__ == '__main__':
    """
    node1 = TreeNode('a')
    node2 = TreeNode('b')
    node3 = TreeNode('c')
    node4 = TreeNode('d')
    node5 = TreeNode('e')
    node6 = TreeNode('f')
    node7 = TreeNode('g')

    node1.left = node2
    node1.right = node3
    node2.left = node4
    node2.right = node5
    node4.left = node6
    node3.left = node7
    """


    list1 = ['a', 'b', 'c', 'd']
    list2 = ['e', 'f', 'g', 'h']

    nodes = {}
    # store the node information
    level = len(list1)
    for i in range(level):
        if i == 0:
            nodes[0] = TreeNode("0")
        else:
            this_level_max = 2 ** (i + 1) - 2
            this_level_min = 2 ** i - 1
            #print("{}\t{}".format(this_level_min, this_level_max))
            for j in range(this_level_min, this_level_max+1):
                if is_odd(j):
                    nodes[j] = TreeNode(list2[i])
                else:
                    nodes[j] = TreeNode(list1[i])

    # store node relationship
    """
                    node1
                    /   \
          node3(left)   node4(right)
    """
    for i in range(level-1):
        if i == 0:
            nodes[0] = TreeNode("0")
            nodes[0].left = nodes[1]
            nodes[0].right = nodes[2]
        else:
            this_level_max = 2 ** (i + 1) - 2
            this_level_min = 2 ** i - 1
            for j in range(this_level_min, this_level_max+1):
                nodes[j].left = nodes[2*j+1]
                nodes[j].right = nodes[2*j+2]


    r = []
    dfs(nodes[0], result=r)
    print(r)

