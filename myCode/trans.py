import os

f = open("./graphviz/twitter.txt")
f2 = open("./graphviz/transed_twitter.txt",'w')
line = f.readline()
num = 0

dic=dict()

while line:
  a = line.split()[0]
  b = line.split()[1]
  if a not in dic:
    dic[a]=num
    num+=1
  if b not in dic:
    dic[b]=num
    num+=1
  print(num)
  f2.writelines(str(dic.get(a))+" "+str(dic.get(b))+"\n")
  line = f.readline()
  
  

f. close()
f2.close()
