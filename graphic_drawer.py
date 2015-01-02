import dnadigest
import pprint
import svgwrite
import math
#After Autosave


def drawer(total_nucl, enzyme_names, nucl_cut_ind):
    fasta_id = '>FASTA #####'
    radius = 250.
    distance = 25
    center = 350
    i=1
    FASTA_headers = ['>FASTA 12345']
    for triad in zip(FASTA_headers, nucl_cut_ind, enzyme_names):
        j=0
        header = triad[0]
        nucl_cut_index = triad[1]
        enzymes = triad[2]
        #total_nucl,  nucl_cut_index, enzyme_names,  and FASTA_headers need to be obtained from the dnadigest class.
        filename = 'plasmid'+str(i)+'.svg'
        i+=1
        svg_document = svgwrite.Drawing(filename = filename, size = ("1000px", "1000px"))
        svg_document.add(svg_document.circle(center = (center, center), r=radius, fill="rgb(255, 255, 255)", stroke="black"))
        svg_document.add(svg_document.text(header, insert = (300,  350)))

        print nucl_cut_index
        for index in nucl_cut_index:
            angle = (float(index)/float(total_nucl))*2*(math.pi)
            line_x = radius*math.cos(angle)+center
            line_y = radius*math.sin(angle)+center
            if index == 0 or index == 0.5*total_nucl:
                slope = 0
            elif index == .25*total_nucl or index == .75*total_nucl:
                slope = 'undefined'
            elif line_y>=center:
                slope = -1*(0.5*(radius**2-line_x**2+2*line_x*center-center**2)**(-0.5)*(2*center-2*line_x))**(-1)
            else:
                slope = -1*(0.5*(radius**2-line_x**2+2*line_x*center-center**2)**(-0.5)*(2*line_x-2*center))**(-1)
            if slope!= 'undefined':
                x1, y1, x2, y2 = line_endpoint_calculator(distance, line_x, line_y, slope)
            else:
                if index == .25*total_nucl:
                    x1 = center
                    x2 = x1
                    y1 = center + radius + distance
                    y2 = center + radius - distance
                else:
                    x1 = center
                    x2 = x1
                    y1 = center - radius + distance
                    y2 = center - radius - distance
            svg_document.add(svg_document.line(start = (x1, y1), end=(x2, y2), stroke="black"))
            text_dist = 25
            if slope == 'undefined':
                if index == .25*total_nucl:
                    x3 = x1
                    y3 = y2 + text_dist
                    svg_document.add(svg_document.text(enzymes[j], insert=(x3, y3)))
                else:
                    x3 = x1
                    y3 = y2 - text_dist
                    svg_document.add(svg_document.text(enzymes[j], insert=(x3, y3)))
            elif line_x >=center and line_y>=center:
                x3, y3, x4, y4 = line_endpoint_calculator(text_dist, x1, y1, slope)
                svg_document.add(svg_document.text(enzymes[j], insert=(x3, y3)))
            elif line_x<=center and line_y>=center:
                x3, y3, x4, y4 = line_endpoint_calculator(text_dist, x2, y2, slope)
                svg_document.add(svg_document.text(enzymes[j], insert=(x4, y4)))
            elif line_x<=center and line_y<=center:
                x3, y3, x4, y4 = line_endpoint_calculator(text_dist, x2, y2, slope)
                svg_document.add(svg_document.text(enzymes[j], insert=(x4, y4)))
            else:
                x3, y3, x4, y4 = line_endpoint_calculator(text_dist, x1, y1, slope)
                svg_document.add(svg_document.text(enzymes[j], insert=(x3, y3)))
            j+=1
        svg_document.add(svg_document.rect(insert = (800,  50), size = ("200px",  "500px"), stroke_width = "1", stroke = "black", fill = "rgb(255, 255, 255)"))
        print(svg_document.tostring())
        svg_document.save()


def line_endpoint_calculator(d, p, q, m):
    x1 = float(p+math.sqrt((d**2)/(m**2+1)))
    y1 = float(m*x1-m*p+q)
    x2 = float(p-math.sqrt((d**2)/(m**2+1)))
    y2 = float(m*x2-m*p+q)
    print 'Point 1: (%s,  %s)' % (x1,  y1)
    print 'Point 2: (%s,  %s)' % (x2,  y2)
    print 'Distance:' + str(math.sqrt((x1-x2)**2+(y1-y2)**2))
    return x1, y1, x2, y2


if __name__ == '__main__':
    digester = dnadigest.Dnadigest()
    enzyme_dict = digester.get_dict('enzyme_data.yaml')
    fragment_list, assoc_enzyme_list, line_marker_list, length = digester.process_data(['AAAAATGTACAAATGTACAAAA'], enzyme_dict, ['AaaI'])
    print 'Fragments:', fragment_list
    print 'Enzymes:', assoc_enzyme_list
    print 'Cut Points:', line_marker_list
    print 'Total Num of Nucleotides:', length
    drawer(length, assoc_enzyme_list, line_marker_list)
