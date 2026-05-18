---
name: baby-teach
description: Use when user need to learn a project(science theory, code, etc.).
---

## Overview

three teaching way for "expert", "general", "baby"

## When to Use

**Explicit triggers:** user says "teach", "do not understand", "不理解", "教教我", "我是笨蛋", "保姆级教学"...

**Proactive triggers:**
- Explain the archtecture and motivation of current project
- How does the equations derative? show more detail

## Workflow

1. ask user: which files needed to be teach? "all"/"special one"/"special some"
2. ask user: mode choose: "expert"(a macro view of current question), "general", "baby"(a detailed answer, the aim is to make rookie user learn enough to rebuild the project and theory by themself)
3. ask user: which language and programming languages you want to use for comparison？(default: english+Chinese｜cpp,python,rust,fortran)
4. ask user: how detailed need to be?(grammer/cn-en noun comparing/compile detail/math foundation...)
5. 将问题拆解成极小的一个个部分（如代码，拆成一个个语法完整的逻辑行）并依次进行讲解，讲解完一行后给出选择等待用户回应（understand/not understand:...）。如问题代码中每个function example 拆成一行一行进行讲解。
默认baby mode情况下要给出专业名词的多语言对比、代码功能在不同程序语言的实现方式对比、实现所用到的语法细节以及数据结构，推导过程中要经常性询问：用户你了解这个数学含义吗？
6. 持续上述流程直到用户对某个部分完全没有问题（即function example这个有完整输入输出的部分或一个推导的逻辑自洽的中间结果）后询问用户：对于...是如何理解的？判断用户答案是否正确，正确则继续接下来的部分，不正确则继续讲解。
7. 回答问题的过程中生成baby-QA-date.md 记录所有问答
