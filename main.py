# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 15:13:23 2019

@author: KurenkovAA
"""
import numpy as np

output = open("out.txt","w",encoding='utf-8')
inputFile = open("in.txt","r")

class Graph:
    
    def addEdge(self,l,r):
        l-=1
        r-=1
        self.out_edge[l] += 1
        self.in_edge[r] += 1
        self.edgesNum[(l,r)] = 2*len(self.edges)
        self.edgesNum[(r,l)] = 2*len(self.edges) + 1
        self.edges.append([l,r])
        self.a[l].append(r)
        self.a[r].append(l)
        self.table = []
        
    def addEdgeToTree(self,e):
        self.treeEdges.append(e)
        [l,r] = self.edges[e]
        self.tree[l].append(r)
        self.tree[r].append(l)
        
    def addEdgeToCotree(self,e):
        self.cotreeEdges.append(e)

    def deleteLastTreeEdge(self):
        e = self.treeEdges[-1]
        [l,r] = self.edges[e]
        self.tree[l].pop()
        self.tree[r].pop()
        self.treeEdges.pop()
        
    def __init__(self):
        
        self.n,self.m = map(int,inputFile.readline().split())
        self.in_edge = [0 for _ in range(self.n)]
        self.out_edge = [0 for _ in range(self.n)]
        self.a = [[] for i in range(self.n + 1)]
        self.tree = [[] for i in range(self.n + 1)]
        self.treeEdges = []
        self.cotreeEdges = []
        self.edgesNum = {}
        self.edges = []
        
        for i in range(self.m):
            self.addEdge(*map(int,inputFile.readline().split()))
            
        ed = [*map(int,inputFile.readline().split())]
        for i in ed:
            self.addEdgeToTree(i-1)
            
        ed = [*map(int,inputFile.readline().split())]
        for i in ed:
            self.addEdgeToCotree(i-1)
        
    def cycleDfs(self,v,p,start):
        for to in self.tree[v]:
            if p == to:
                continue
            if to == start:
                return [v]
            b = self.cycleDfs(to,v,start)
            if len(b) != 0:
                b.append(v)
                return b
        return []
    
    def getCycle(self,e,index):
        [l,_] = self.edges[e]
        self.addEdgeToTree(e)
        
        cycle = self.cycleDfs(l,-1,l)
        vector = [0 for i in range(self.m)]
        self.deleteLastTreeEdge()
        
        output.write(u"c_{}⇒ℂ_{}=".format(index + 1,index + 1))
        output.write("{")
        
        for i in range(len(cycle)):
            res = self.edgesNum[cycle[i-1],cycle[i]]
            num = res//2
            sign = (1 if res % 2 == 0 else -1)
            vector[num] += sign
            if i != 0:
                output.write(",")
            output.write("a_{}".format(num+1))
        output.write("}")
        output.write("~I_{}=(".format(index + 1))
        
        self.table.append(vector)
        for i in range(len(vector)):
            if i != 0:
                output.write(",")
            output.write(str(vector[i]))
            
        output.write(")^T\n")
        
    def getCycles(self):
        for i,val in enumerate(self.cotreeEdges):
            self.getCycle(val,i)
            
            
    def cocycleDfs(self,v,p,ed):
        res = [0 for i in range(self.m)]
        
        for to in self.a[v]:
            e = self.edgesNum[(to,v)]
            num = e//2
            sign = (1 if e % 2 == 0 else -1)
            res[num] += sign
        
        for to in self.tree[v]:
            e = self.edgesNum[(to,v)]
            num = e//2          
            if p == to or num == ed:
                continue
            b = self.cocycleDfs(to,v,ed)
            res = [x + y for (x,y) in zip(res,b)]
            
        return res
            
    
    def getCocycle(self,e,index):
        [_,r] = self.edges[e]
        
        vector = self.cocycleDfs(r,-1,e)
        self.table.append(vector)
        
        output.write(u"b_{}⇒S_{}=".format(index + 1,index + 1))
        output.write("{")
        
        k = 0
        for i,val in enumerate(vector):
            if val == 0:
                continue
            if k != 0:
                output.write(",")
            k += 1
            output.write("a_{}".format(i+1))
        output.write("}")
        output.write("~U_{}=(".format(index + 1))
        
        for i in range(len(vector)):
            if i != 0:
                output.write(",")
            output.write(str(vector[i]))
            
        output.write(")^T\n")
            
    def getCocycles(self):
        for i,val in enumerate(self.treeEdges):
            self.getCocycle(val,i)
            
    def out(self):
        self.getCycles()
        output.write("\n")
        self.getCocycles()
        output.write("\n")
        for row in self.table:
            for val in row:
                output.write(u"{}\t".format(val))
            output.write("\n")
    
    def check(self):
        for ind,(in_num,out_num) in enumerate(zip(self.in_edge,self.out_edge)):
            print(ind + 1,in_num,out_num)

g = Graph()
g.out()

r = [float(val) for val in inputFile.readline().split()]
print(r)
E = np.array([float(val) for val in inputFile.readline().split()])
print(E)
M = np.transpose(np.array(g.table))

# r = np.diag([1,2])
# e = np.array([3,6])
# M = r @ np.array([[2,1],[1,2]])

v,e = g.n,g.m
L = np.array([[res*val if ind < e-v + 1 else -val
      for (ind,val) in enumerate(line)] for (line,res) in zip(M,r)])
print(M)
print(L)

# print(r)
# print(e)

Mi = M.T[:e - v + 1,:]
Mu = M.T[e - v + 1:,:]


g.check()
x = np.linalg.solve(L,E.T)
I = x[:e - v + 1] @ Mi
U = x[e - v + 1:] @ Mu

def monomial(coeff,ind,line):
    coeff = int(coeff)
    ind = int(ind)
    fl = sum([val != 0 for val in line[:ind-1]]) == 0
    if (coeff == 0):
        return ''
    if (coeff == 1):
        return ('x_{}' if fl else '+x_{}').format(ind)
    if (coeff == -1):
        return '-x_{}'.format(ind)
    return ('{}x_{}' if fl or coeff < 0 else '+{}x_{}').format(coeff,ind)

def sysEqString(L,r):
    return '{█(' + '@'.join(['&'.join([monomial(val,ind + 1,line) for ind,val in enumerate(line)]) + '&=&{}'.format(int(eq)) for line,eq in zip(L,r)]) + ')┤'

output.write(sysEqString(L,E) + '\n')
output.write('X=(■({})) '.format('@'.join([f'{val:.2f}' for val in x])))
output.write('I=(■({})) '.format('@'.join([f'{val:.2f}' for val in I])))
output.write('U=(■({}))\n'.format('@'.join([f'{val:.2f}' for val in U])))
inputFile.close()
output.close()
