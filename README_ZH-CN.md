nsapy
=====

A Python Module for Nonlinear Structural Analysis

编号规则
-------
	* 每种对象都有标签(tag)属性，多个对象实例组成字典(Dict)时，标签(Tag)属性作为对象字典(Dict)的检索键(Keys)，值(Values)为对象实例本身。
	* 节点具有索引(index)属性，用于定位其在整体刚度、质量矩阵中的位置。节点的索引(index)值按节点的输入顺序依次定义为 0,1,2,3...
	* 所有节点输入完成后，根据节点的index,dof值确定待求解的自由度编号序列，此序列为隐含变量，可表示为range(sum([node_i.dof]))