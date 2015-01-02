import svgwrite
import math
def drawer():
    fasta_id = '>FASTA #####'
    radius = 250.
    total_nucl = 10000.
    nucl_cut_index = 0.
    distance = 500
    angle = (nucl_cut_index/total_nucl)*2*(math.pi)
    line_x = radius*math.cos(angle)
    line_y = radius*math.sin(angle)
    if nucl_cut_index == 0 or nucl_cut_index == 0.5*total_nucl:
        slope = 0
    elif (line_x>=0 and line_y>=0) or (line_x<=0 and line_y>=0):
        slope = -1*1/(-1*radius*math.cos(angle)*math.sqrt(radius**2-(radius**2)*(math.cos(angle))**2))
    else:
        slope = -1*1/(radius*math.cos(angle)*math.sqrt(radius**2-(radius**2)*(math.cos(angle))**2))

    svg_document = svgwrite.Drawing(filename = "test-svgwrite.svg",size = ("1000px","1000px"))
    svg_document.add(svg_document.circle(center = (350,350),r=radius,fill = "rgb(255,255,255)",stroke="black"))
    svg_document.add(svg_document.line(start = (line_x+distance,slope*(line_x+distance)-slope*line_x-line_y),end = (line_x-distance,slope*(line_x-distance)-slope*line_x-line_y)))
    print line_x+distance,slope*(line_x+distance)-slope*line_x-line_y
    print line_x-distance,slope*(line_x-distance)-slope*line_x-line_y
    svg_document.add(svg_document.text(fasta_id,insert = (300, 350)))
    print(svg_document.tostring())
    svg_document.save()

def line_endpoint_calculator():
    d = float(2)
    p = float(0)
    q = float(0)
    m = float(math.sqrt(3))
    x1 = float(p+math.sqrt((d**2+(m*p)**2+(p*m**2)**2)/(m**2+1)))
    y1 = float(m*x1-m*p+q)
    x2 = float(p-math.sqrt((d**2+(m*p)**2+(p*m**2)**2)/(m**2+1)))
    y2 = float(m*x2-m*p+q)
    print 'Point 1: ('+str(x1)+','+str(y1)+')'
    print 'Point 2: ('+str(x2)+','+str(y2)+')'
    print 'Twice given distance:',2.0*d
    print 'Calculated distance 1:',math.sqrt((x1-p)**2+(y1-q)**2)
    print 'Calculated distance 2:',math.sqrt((x2-p)**2+(y2-q)**2)



