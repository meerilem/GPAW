import math
#from gi.repository import cairo
import gizeh
from matplotlib import rc

ncolors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
ncolors_norm=[(0.216, 0.494, 0.722), (1, 0.498, 0), (0.302, 0.686, 0.29), (0.969, 0.506, 0.749), (0.651, 0.337, 0.157), (0.596, 0.306, 0.639)]

rc('text', usetex=True)
Pi = 3.14
ls=120
i_s=80
center_x=500+10
max_y=570
ff="Impact"
rx=10
ry=200
surface = gizeh.Surface(width=1310, height=max_y) # in pixels
rect1 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x,max_y-60), fill=(0,0,0), angle=Pi/2)
rect2 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x,max_y-40), fill=(0,0,0), angle=Pi/2)
text1 = gizeh.text(r'C', fontfamily=ff,  fontsize=ls, fill=ncolors_norm[1], xy=(center_x-170, max_y-50))
text1_2 = gizeh.text(r'2', fontfamily=ff,  fontsize=i_s, fill=ncolors_norm[1], xy=(center_x-120, max_y-30))
text2 = gizeh.text(r"C", fontfamily=ff,  fontsize=ls, fill=ncolors_norm[2], xy=(center_x+150, max_y-50))
text2_2 = gizeh.text(r"3", fontfamily=ff,  fontsize=i_s, fill=ncolors_norm[2], xy=(center_x+200, max_y-30))
rect3 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x-200,max_y-195), fill=(0,0,0), angle=(-3*Pi/5-Pi/2))
rect4 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x+200,max_y-195), fill=(0,0,0), angle=(3*Pi/5-Pi/2))
text3 = gizeh.text("N", fontfamily=ff,  fontsize=ls, fill=(0,0,0), xy=(center_x-240,max_y-340))
text4 = gizeh.text("N", fontfamily=ff,  fontsize=ls, fill=(0,0,0), xy=(center_x+240,max_y-340))
text3_2 = gizeh.text("1", fontfamily=ff,  fontsize=i_s, fill=(0,0,0), xy=(center_x-185,max_y-310))
text4_2 = gizeh.text("2", fontfamily=ff,  fontsize=i_s, fill=(0,0,0), xy=(center_x+295,max_y-310))
rect5 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x-120,max_y-440), fill=(0,0,0), angle=(-6*Pi/5-Pi/2))
rect6 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x+140,max_y-447), fill=(0,0,0), angle=(6*Pi/5-Pi/2))
rect7 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x-134,max_y-454), fill=(0,0,0), angle=(-6*Pi/5-Pi/2))
text5 = gizeh.text(r'C', fontfamily=ff,  fontsize=ls, fill=ncolors_norm[0], xy=(center_x,max_y-520))
text5_2 = gizeh.text(r'1', fontfamily=ff,  fontsize=i_s, fill=ncolors_norm[0], xy=(center_x+50,max_y-500))
rect8 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x-360,max_y-414), fill=(0,0,0), angle=(6*Pi/5-Pi/2))
rect9 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x+360,max_y-414), fill=(0,0,0), angle=(-6*Pi/5-Pi/2))
text6 = gizeh.text(r'C', fontfamily=ff,  fontsize=ls, fill=ncolors_norm[3], xy=(center_x-470,max_y-520))
text6_2 = gizeh.text(r'4', fontfamily=ff,  fontsize=i_s, fill=ncolors_norm[3], xy=(center_x-410,max_y-500))
text7 = gizeh.text(r"C", fontfamily=ff,  fontsize=ls, fill=ncolors_norm[4], xy=(center_x+455,max_y-520))
text7_2 = gizeh.text(r"5", fontfamily=ff,  fontsize=i_s, fill=ncolors_norm[4], xy=(center_x+505,max_y-500))
rect10 = gizeh.rectangle(lx=rx, ly=ry, xy=(center_x+600,max_y-414), fill=(0,0,0), angle=(6*Pi/5-Pi/2))
text8 = gizeh.text("C", fontfamily=ff,  fontsize=ls, fill=ncolors_norm[5], xy=(center_x+720,max_y-340))
text8_2 = gizeh.text("6", fontfamily=ff,  fontsize=i_s, fill=ncolors_norm[5], xy=(center_x+780,max_y-310))
text9 = gizeh.text("+", fontfamily=ff,  fontsize=i_s, fill=(0,0,0), xy=(center_x-240,max_y-410))

rect1.draw(surface) 
rect2.draw(surface) 
text1.draw(surface)
text1_2.draw(surface)
text2.draw(surface)
text2_2.draw(surface)
rect3.draw(surface)
rect4.draw(surface)
text3.draw(surface)
text4.draw(surface)
text3_2.draw(surface)
text4_2.draw(surface)
rect5.draw(surface)
rect6.draw(surface)
rect7.draw(surface)
text5.draw(surface)
text5_2.draw(surface)
rect8.draw(surface)
rect9.draw(surface)
text6.draw(surface)
text7.draw(surface)
text6_2.draw(surface)
text7_2.draw(surface)
rect10.draw(surface)
text8.draw(surface)
text8_2.draw(surface)
text9.draw(surface)
surface.write_to_png("EMIm.png") # export the surface as a PNG
