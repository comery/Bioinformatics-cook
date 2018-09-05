import os
import brewer2mpl
import argparse
import subprocess
import string

parser = argparse.ArgumentParser(
	description = "Draw circos figure to illustrat mapping result")

parser.add_argument('-fa', metavar = 'FILE', type = str, required = True,
                    dest = 'ref', help = "reference fasta")
parser.add_argument('-map', metavar = 'FILE', type = str, required = True,
		    help = "mapping result, depth of each site")
parser.add_argument('-conf', metavar = 'FILE', type = str, required = True,
		    help = "general configures to circos parameters")

args = parser.parse_args()

try:
    fh=open(args.ref)
except Exception as e:
    print("can not find file " + args.ref)
    exit(0)

''' read configures '''

configures = {}
with open(args.conf, 'r') as fh_conf:
    for c in fh_conf.readlines():
        c = c.strip()
        if len(c) >0 and c[0] != "#":
            vals = c.split("=")
            configures[vals[0].strip()] = vals[1].strip()

print("read configures file done!")
if os.path.exists(configures['outdir']) == False:
    os.mkdir(configures['outdir'])

seq = {}
chrs_id = {}
for line in fh:
    if line.startswith('>'):
        name=line.replace('>','')
        name= name.strip()
        seq[name] = ""
    else:
        seq[name]+= line.replace('\n','')
fh.close()

num = 0
for k,v in seq.items():
    seq[k] = len(v) # change original id , sequence pairs into id-lenght
    num += 1
    chrs_id[k] = num

fh.close()

chrs_in_depth = []

''' depth feature '''

max_depth = 0
fh_depth = open(configures['outdir'] + "/depth.txt", 'w')
try:
    fh_map = open(args.map,'r')
except Exception as e:
    print(e)
    exit(0)
for m in fh_map:
    a = m.split("\t")
    name = a[0]
    if name not in chrs_in_depth:
        chrs_in_depth.append(name)
    pos = str(a[1])
    dep = a[2]
    if int(dep) > max_depth:
        max_depth = int(dep)
    fh_depth.write("hs" + str(chrs_id[name])  + " " + pos + " " + pos + " " + dep )

fh_depth.close()
refs_num = len(chrs_in_depth)

''' colors theme '''
colors = []
if refs_num < 3:
    bmap = brewer2mpl.get_map('Paired', 'Qualitative', 3, reverse=True)
    colors = bmap.colors

elif refs_num >= 3 and refs_num < 13:
    bmap = brewer2mpl.get_map('Paired', 'Qualitative', 12, reverse=True)
    colors = bmap.colors
elif refs_num >= 13 and refs_num < 34:
    brbg = brewer2mpl.get_map('BrBG', 'Diverging', 11, reverse=True)
    colors = brbg.colors
    piyg = brewer2mpl.get_map('PiYG', 'Diverging', 11, reverse=True)
    colors += piyg.colors
    paired = brewer2mpl.get_map('Paired', 'Qualitative', 12, reverse=True)
    colors += paired.colors
else:
    print("Sorry, too many chrs, I can not deal that!")
    exit(0)

'''print("colors for usage: " + colors)'''

''' karyotype '''

fh_karyo = open(configures['outdir'] + '/karyotype.txt', 'w')

for k in range(0,refs_num):
    this_color = ','.join(str(v) for v in colors[k])
    fh_karyo.write("chr" + str(k+1) + " - " +  "hs" + str(k+1) + " " + chrs_in_depth[k] +\
                    " 0 " + str(seq[chrs_in_depth[k]]) + " " + this_color + "\n")
fh_karyo.close()


''' generate final circos configure '''

with open(configures['outdir'] + "/circos.conf",'w') as circos:
    circos.write('''
<<include etc/colors_fonts_patterns.conf>>

#-----------------image------------------
<image>
###<<include etc/image.conf>>
''')
    circos.write("dir   =  " + configures['outdir'] + "\n")
    circos.write("file  =  circos.png" +"\n")
    circos.write("png   =  " + configures['png'] + "\n")
    circos.write("svg   =  " + configures['svg'] + "\n")
    circos.write('''
# radius of inscribed circle in image
radius         = 1500p

# by default angle=0 is at 3 o'clock position
angle_offset      = -90

#angle_orientation = counterclockwise
auto_alpha_colors = yes
auto_alpha_steps  = 5
''')
    circos.write("background = " + configures['background'] + "\n")
    circos.write('''
</image>

#-----------------ideogram------------------
<ideogram>

<spacing>
default = 0.01r
break   = 0.1r
</spacing>

###<<include ideogram.position.conf>>
radius           = 0.80r
thickness        = 20p
fill             = yes
fill_color       = grey
stroke_thickness = 3
stroke_color     = black

###<<include ideogram.label.conf>>
show_label       = yes
label_font       = bolditalic
label_radius     = dims(image,radius) + 0.01r
#label_radius     = 0.8r + 40p
label_size       = 15p
label_parallel   = yes
label_case       = lowwer
#label_format     = eval(sprintf("chr%s",var(label)))
#label_format     = eval(var(labe))


###<<include bands.conf>>
show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 0

</ideogram>
#-----------------ticks------------------
show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = dims(ideogram,radius_outer)
#radius           = 1r+0.06r
orientation      = out
label_multiplier = 1e-3
color            = black
thickness        = 2p
font             = blod

<tick>
spacing        = 1u
show_label     = no
label_size     = 15p
size           = 15p
format         = %d
label_offset   = 2p
#suffix         = " kb"
</tick>

<tick>
spacing        = 5u
show_label     = yes
label_size     = 20p
size           = 20p
format         = %d
suffix         = " kb"
label_offset   = 2p
</tick>

<tick>
spacing        = 10u
show_label     = yes
label_size     = 20p
size           = 20p
format         = %d
label_offset   = 2p
suffix         = " kb"
</tick>

</ticks>
#-----------------karyotype------------------
''')
    circos.write("karyotype   = " + configures['outdir'] + "/karyotype.txt" + "\n")

    circos.write('''
chromosomes_units = 1000
chromosomes_display_default = yes

#-----------------plots------------------

<plots>

########### depth information
<plot>
type      = line
thickness = 2
max_gap = 1u
''')

    circos.write("file    = " + configures['outdir'] + "/depth.txt" + "\n")
    circos.write("color   = dgreen" + "\n")
    circos.write("min     = 0" + "\n")
    circos.write("max     = " + str(max_depth) + "\n")
    circos.write("r0      = 0.5r" + "\n" )
    circos.write("r1      = 0.68r" + "\n")
    circos.write("fill_color = " + configures['depth_fill'] + "\n")

    circos.write('''
<axes>
<axis>
color     = lgrey_a2
thickness = 1
spacing   = 0.06r
</axis>
</axes>

<rules>

<rule>
condition    = var(value) > 5
color        = dred
fill_color   = dred_a1
</rule>

#<rule>
#condition    = var(value) < 3
#color        = dred
#fill_color   = dred_a1
#</rule>

</rules>

</plot>
</plots>

<<include etc/housekeeping.conf>>

'''
)

    print("Drawing circos.\n")
cmd = configures['circos_path'] + " -conf " + configures['outdir'] + "/circos.conf"
subprocess.call(cmd, shell = True)

print("All done !")
