try:
  import itchat,pyecharts
except ImportError:
  mag = "you did not install module - itchat, and also you need to check other module:\n\t" +\

        "matplotlib.pyplot\n" +\
        "pyecharts\n" +\
        "echarts-countries-pypkg\n" +\
        "echarts-china-provinces-pypkg\n" +\
        "echarts-china-cities-pypkg\n" +\
        "echarts-china-counties-pypkg\n" +\
        "echarts-china-misc-pypkg\n" 


import itchat               # 登陆微信，获取朋友信息的库
import collections           # 计数用的库
import matplotlib.pyplot as plt # 画饼图用的库
from matplotlib.font_manager import FontProperties # 解决图片中文乱码
from pyecharts import Map       # 画地图用的库
from pyecharts import Geo       # 画地图二

'''
pip install echarts-countries-pypkg    #全球国家地图: echarts-countries-pypkg (1.9MB): 世界地图和 213 个国家，包括中国地图
pip install echarts-china-provinces-pypkg # 中国省级地图: echarts-china-provinces-pypkg (730KB)：23 个省，5 个自治区
pip install echarts-china-cities-pypkg  # 中国市级地图: echarts-china-cities-pypkg (3.8MB)：370 个中国城市
pip install echarts-china-counties-pypkg # 中国县区级地图: echarts-china-counties-pypkg (4.1MB)：2882 个中国县·区
pip install echarts-china-misc-pypkg   # 中国区域地图: echarts-china-misc-pypkg (148KB)：11 个中国区域地图，比如华南、华北
'''

itchat.login()          # 登陆微信网页版
friends = itchat.get_friends(update=True)[0:] # 获取当前微信好友信息
print(friends)
sexs = list(map(lambda x: x['Sex'],friends[1:]))

# 解释下这行代码，微信中性别字段的取值有Unkonw、Male和Female三种，其对应的数值分别为0、1、2。
# 我们将会得到如 sexs = [1,0,1,0,1,0,1,0,2,1,1,2,1] 的集合
# 这行代码语法看不懂的，可以对应看这下面行代码帮助语法上的理解 ，
# map(lambda x: x ** 2, [1, 2, 3, 4, 5])  会得到 [1, 4, 9, 16, 25]

sex_counts = [0,1,2]
sex = collections.Counter(sexs) # 通过Collection模块中的Counter()对这三种不同的取值进行统计
counts = []                     # 性别统计结果
for i in sex_counts:            # 按照 0，1，2 的顺序统计出相应的性别，
   counts.append(sex[i])
print(counts)
labels = ['Unknow','Male','Female'] # 设置饼图的标签
colors = ['red','blue','coral']     # 设置饼图的颜色
explode = [0,0,0.1]                 # 0.1 凸出这部分，female
plt.figure(figsize=(8, 8))          # 设置绘图对象的大小
plt.axes(aspect=1)                  # 参数设置这个饼为正圆
plt.pie(counts,                     # 性别统计结果
       labels=labels,              # 性别展示标签
       colors=colors,              # 饼图区域配色
       explode=explode,            # 凸出部分
       labeldistance=1.1,          # 标签距离圆点距离，1.1指1.1倍半径的位置
       autopct='%3.1f%%',          # 饼图区域文本格式,%3.1f%%表示小数有三位，整数有一位的浮点数
       shadow=False,               # 饼图是否显示阴影
       startangle=90,              # 起始角度，0，表示从0开始逆时针转，为第一块。一般选择从90度开始比较好看
       pctdistance=0.6             # 饼图区域文本距离圆点距离
       )
font_set = FontProperties(fname=r"c:\windows\fonts\simsun.ttc", size=15) # 解决 Windows 环境下乱码问题
plt.title(u'Fightjiang 的微信好友性别比例', fontproperties=font_set)
plt.show()

attr = ['安徽', '北京', '福建', '广东', '贵州', '海南', '河北', '河南', '黑龙江',
       '湖北', '湖南', '吉林', '江苏', '辽宁', '山东','山西', '陕西', '上海',
       '四川', '天津', '云南', '浙江', '重庆']
friend = []    # 好友所在的省份
for i in friends[1:]:
   friend.append(i['Province'])
print(friend)
location = collections.Counter(friend) #  一个迭代对象生成的counter

value = []     # weixin1 Map
for i in attr:  # value 每个省会对应的数量
   value.append(location[i])
print(value)
map = Map(itchat.search_friends()['NickName'] +  " 各省微信好友分布",width=1200, height=600)
map.add("", attr, value, maptype='china', is_visualmap=True,
       visual_text_color='#000')
#map.show_config()
map.render('weixin1.html')

City = [] # 微信好友所在城市
for city in friends[1:]:
   City.append(city['City'])
Citys = collections.Counter(City) # 每个城市对应的数量

values = []                  # weixin2 Map
for city in set(City):           # values 每个城市对应的数量
    if(city != '' and city.isalpha()): # 除去没有城市的 和 外国城市
       values.append((city,Citys[city]))
print(values)

geo = Geo(itchat.search_friends()['NickName'] + " 各省微信好友分布",
         title_color="#fff", title_pos="center",
         width=1200, height=600, background_color='#404a59')
attr, value = geo.cast(values)
geo.add("", attr, value, visual_range=[0, 200],
       visual_text_color="#fff", symbol_size=15, is_visualmap=True)
#  geo.show_config()
geo.render("weixin2.html")


