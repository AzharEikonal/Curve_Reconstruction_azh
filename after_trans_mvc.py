import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point, Polygon
from shapely.geometry.base import CoordinateSequence
from shapely.ops import polygonize, unary_union
import numpy as np
import random
import os
import sys

def parse_obj(file_path):
    vertices = []
    faces=[]
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('v '):
                parts = line.split()
                vertices.append((float(parts[1]), float(parts[2])))
            elif line.startswith('f '):
                parts = line.split()
                face = [int(part.split('/')[0]) -1 for part in parts[1:]]
                faces.append(face)
    return vertices, faces

def build_lines(vertices, edges):
    lines = []
    for edge in edges:
        line = LineString([vertices[edge[0]], vertices[edge[1]]])
        lines.append(line)
    return lines

def plot_polygons(polygons, mesh_polygons):
    fig, ax = plt.subplots()
    colors = plt.cm.get_cmap('hsv', len(polygons))

    for i, poly in enumerate(polygons):
        #poly is a list of coordinates
        x, y = zip(*poly)
        ax.plot(x, y, color=colors(i))
        
        
    for i, mesh_points in enumerate(mesh_polygons):
        x, y = zip(*mesh_points)
        ax.scatter(x, y, color=colors(i), label=f'Polygon {i+1}')
    
    plt.legend()
    plt.show()

def find_intersection_points(rectangle, cutting_lines):
    # Define rectangle boundary lines
    rect_lines = [
        LineString([rectangle[0], rectangle[1]]),
        LineString([rectangle[1], rectangle[2]]),
        LineString([rectangle[2], rectangle[3]]),
        LineString([rectangle[3], rectangle[0]])
    ]

    # Combine all lines (rectangle boundaries + cutting lines)
    all_lines = rect_lines + cutting_lines

    # Find intersection points
    points = []
    for line1 in all_lines:
        for line2 in all_lines:
            if line1 != line2:
                if line1.intersects(line2):
                    intersect = line1.intersection(line2)
                    if intersect.geom_type == "Point":
                        points.append(intersect)
                    elif intersect.geom_type == "MultiPoint":
                        points.extend(list(intersect))
    
    # Remove duplicates and sort the points
    unique_points = list(set(points))
    return unique_points, all_lines

def form_polygons(intersection_points, lines):
    # Add intersection points as vertices to lines
    new_lines = []
    for line in lines:
        new_coords = [line.coords[0]]
        for point in intersection_points:
            if line.distance(point) < 1e-8 and point.coords[0] not in new_coords:
                new_coords.append(point.coords[0])
        new_coords.append(line.coords[1])
        new_coords = sorted(new_coords)
        for i in range(len(new_coords) - 1):
            new_lines.append(LineString([new_coords[i], new_coords[i+1]]))

    # Form polygons from the lines
    merged = unary_union(new_lines)
    polygons = list(polygonize(merged))
    return polygons

#return indices 
def extract_mesh_for_polygons(polygons, vertices):
    mesh_polygons = []
    for poly in polygons:
        points_in_poly = [vertex for vertex in vertices if poly.contains(Point(vertex)) or poly.touches(Point(vertex))]
        mesh_polygons.append(points_in_poly)
    return mesh_polygons

def read_sdt_file(sdt_file_path):
    sdt_data = {}
    with open(sdt_file_path, 'r') as file:
        for line in file:
            parts = line.split()
            index = int(parts[0])
            distance = float(parts[1])
            sdt_data[index] = distance
    return sdt_data

if len(sys.argv) < 2:
    print("Please provide a folder name as an argument.")
    sys.exit(1)

folder_name = sys.argv[1]
filename_obj = folder_name + '.obj'
file_path_obj = os.path.join(folder_name, filename_obj)
filepath_bbox_poly = os.path.join(folder_name, 'bbox_poly.txt')

rectangle = []
with open(filepath_bbox_poly, 'r') as file:
    for line in file:
        parts = line.split()
        rectangle.append((float(parts[0]), float(parts[1])))
rectangle = [rectangle[0], rectangle[1], rectangle[2], rectangle[3]]

filepath_cutting_lines = os.path.join(folder_name, 'cutting_lines.txt')
cutting_lines = []
with open(filepath_cutting_lines, 'r') as file:
    for line in file:
        parts = line.split()
        cutting_lines.append(LineString([(float(parts[0]), float(parts[1])), (float(parts[2]), float(parts[3]))]))
                      

# Parse the OBJ file (example file path)
file_path= file_path_obj
vertices, faces = parse_obj(file_path)

# Find intersection points
intersection_points, all_lines = find_intersection_points(rectangle, cutting_lines)

# Form polygons
polygons = form_polygons(intersection_points, all_lines)

# read vertices from sdt file and push them into an array
sdt_name = folder_name + '.sdt'
filepath_sdt = os.path.join(folder_name, sdt_name)
vertices_on_cutting_lines = []
with open(filepath_sdt, 'r') as file:
    for line in file:
        parts = line.split()
        index = int(parts[0])
        vertices_on_cutting_lines.append(vertices[index])

# #plot vertices on cutting lines as scatter plot
# fig, ax = plt.subplots()
# x, y = zip(*vertices_on_cutting_lines)
# ax.scatter(x, y, color='red', label='Vertices on Cutting Lines')

# plt.show()

# Extract mesh for each polygon
mesh_polygons = extract_mesh_for_polygons(polygons, vertices)

# convert polygons into a list of coordinates
polygons = [list(poly.exterior.coords) for poly in polygons]

# calculate minimum distance betweeen vertices of mesj
min_distance = 1000000
for i in range(len(vertices)):
    for j in range(i+1, len(vertices)):
        distance = np.sqrt((vertices[i][0]-vertices[j][0])**2 + (vertices[i][1]-vertices[j][1])**2)
        if distance < min_distance:
            min_distance = distance

# # find vertices from vertices_on_cutting_lines lying on the boundary of first polygon 
# boundary_vertices = []
# for i in range(len(polygons[0])):
#     #make a line from ith vertex to i+1th vertex of the polygon and use distance between line and point. Do not use line string
#     line=[polygons[0][i], polygons[0][(i+1)%len(polygons[0])]]
#     for vertex in vertices_on_cutting_lines:
#         # use line and point distance formula and vertex must lie between the two points of the line
#         # first calculate distance between line and point
#         # then calculate distance between two points of the line
#         # then calculate distance between point and two points of the line
#         # if sum of the two distances is equal to distance between two points of the line then point lies on the line
#         # then check if the point lies between the two points of the line
#         # if it does then add it to the boundary vertices
#         distance_line_point = abs((line[1][1]-line[0][1])*vertex[0] - (line[1][0]-line[0][0])*vertex[1] + line[1][0]*line[0][1] - line[1][1]*line[0][0]) / np.sqrt((line[1][1]-line[0][1])**2 + (line[1][0]-line[0][0])**2)
#         distance_line = np.sqrt((line[1][1]-line[0][1])**2 + (line[1][0]-line[0][0])**2)
#         distance_point_line = np.sqrt((line[0][1]-vertex[1])**2 + (line[0][0]-vertex[0])**2) + np.sqrt((line[1][1]-vertex[1])**2 + (line[1][0]-vertex[0])**2)
#         if distance_line_point < min_distance and abs(distance_line-distance_point_line) < min_distance:
#             boundary_vertices.append(vertex)
        

# #plot vertices on the boundary of first polygon as scatter plot
# fig, ax = plt.subplots()
# x, y = zip(*boundary_vertices)
# ax.scatter(x, y, color='green', label='Vertices on Boundary of First Polygon')
# plt.legend()
# plt.show()

# Now update array of mesh_polygons to include boundary vertices do above for all polygons
for i, poly in enumerate(mesh_polygons):
    boundary_vertices = []
    for j in range(len(polygons[i])):
        #make a line from ith vertex to i+1th vertex of the polygon and use distance between line and point. Do not use line string
        line=[polygons[i][j], polygons[i][(j+1)%len(polygons[i])]]
        for vertex in vertices_on_cutting_lines:
           #use above tactic
            #also put if condition that vertex is not already in the line
            if vertex not in line and vertex not in boundary_vertices and vertex not in poly:
                #calculate distnace between two points of the line
                distance_line = np.sqrt((line[1][1]-line[0][1])**2 + (line[1][0]-line[0][0])**2)
                # calculate two distances of vertex from two points of the line
                d1=np.sqrt((line[0][1]-vertex[1])**2 + (line[0][0]-vertex[0])**2)
                d2=np.sqrt((line[1][1]-vertex[1])**2 + (line[1][0]-vertex[0])**2)
                # if sum of the two distances is equal to distance between two points of the line then point lies on the line
                if abs(distance_line-(d1+d2)) < min_distance:
                    boundary_vertices.append(vertex)    
    mesh_polygons[i] += boundary_vertices


# # plot nth polygonas scattered plot
# n=3
# # plot the nth polygon as scatter plot and plot interior points as red and points on the boundary as green
# fig, ax = plt.subplots()
# x, y = zip(*mesh_polygons[n])
# ax.scatter(x, y, color='red', label='Interior Points')
# plt.legend()
# plt.show()

#now convert mesh_polygon_n array into indices instead of coordinates (indices from original mesh file)
#find indices of vertices on every mesh_polygon_n from vertices

for i, poly in enumerate(mesh_polygons):
    string  = f'map_for_indices_{i} = dict()'
    exec(string)
    for j, vertex in enumerate(poly):
        string = f'map_for_indices_{i}[vertices.index(vertex)] = {j}'
        exec(string)


# calculate indices_mesh_polygon_n for every mesh_polygon_n
for i, poly in enumerate(mesh_polygons):
    string  = f'indices_mesh_polygon_{i} = [' 
    for vertex in poly:
        string += f'{vertices.index(vertex)}, '
    string += ']'
    exec(string)

# # print indices_mesh_polygon_n for every mesh_polygon_n
# for i in range(len(mesh_polygons)):
#     exec(f'print(indices_mesh_polygon_{i})')


# indices_mesh_polygon_n = []
# for i, vertex in enumerate(mesh_polygons[n]):
#     indices_mesh_polygon_n.append(vertices.index(vertex))



# get all those triangles from faces which contains all the indices (if all 3 elements are in indices then add the triangle)
# mesh_polygons_faces = []
# for face in faces:
#     face1=set(face)
#     if face1.issubset(set(indices_mesh_polygon_n)):
#         mesh_polygons_faces.append(face)

# calculate mesh_polygons_faces for every mesh_polygon_n by using intersection of indices_mesh_polygon_n and faces
# find intersection of indices_mesh_polygon_n and faces

for i, poly in enumerate(mesh_polygons):
    string  = f'mesh_polygons_faces_{i} = []'
    exec(string)
    for face in faces:
        face1=set(face)
        if face1.issubset(set(eval(f'indices_mesh_polygon_{i}'))):
            exec(f'mesh_polygons_faces_{i}.append(face)')  


# # print mesh_polygons_faces for every mesh_polygon_n
# for i in range(len(mesh_polygons)):
#     exec(f'print(mesh_polygons_faces_{i})')

# # now find intersection of mesh_polygons_n and vertices_on_cutting_lines 
# boundary_mesh_polygon_n=[]
# for vertex in mesh_polygons[n]:
#     if vertex in vertices_on_cutting_lines:
#         boundary_mesh_polygon_n.append(vertex)

# find boundary_mesh_polygon_n for every mesh_polygon_n
for i, poly in enumerate(mesh_polygons):
    string  = f'boundary_mesh_polygon_{i} = []'
    exec(string)
    for vertex in poly:
        if vertex in vertices_on_cutting_lines:
            exec(f'boundary_mesh_polygon_{i}.append(vertex)')


# make an array distances_mesh_polygon_i for evrry mesh_polygon_i
for i, poly in enumerate(mesh_polygons):
    string = f'distances_mesh_polygon_{i} = []'
    exec(string)
    filepath_mesh_sdt = os.path.join(folder_name, f'distances_mesh_polygon_{i}.txt')
    with open(filepath_mesh_sdt, 'r') as file:
        for line in file:
            values = line.split()
            for value in values:
                if value != 'nan':
                    string = f'distances_mesh_polygon_{i}.append({value})'
                    exec(string)
                else:
                    string = f'distances_mesh_polygon_{i}.append(1000000)'
                    exec(string)
    exec(f'distances_mesh_polygon_{i} = np.array(distances_mesh_polygon_{i})')
    
# update sdt_data dictionary using indices_mesh_polygon_n and distances_mesh_polygon_n
sdt_data= read_sdt_file(filepath_sdt)
for i, poly in enumerate(mesh_polygons):
    for j, val in enumerate(eval(f'distances_mesh_polygon_{i}')):
        #find the values of jth index in indices_mesh_polygon_n
        index= eval(f'indices_mesh_polygon_{i}[j]')
        vertex= vertices[index]
        # check if the vertex is in boundary_mesh_polygon_n
        if vertex not in eval(f'boundary_mesh_polygon_{i}'):
            # update the sdt_data dictionary
            sdt_data[index]=val



# arrange keys of sdt_data in ascending order
sdt_data = dict(sorted(sdt_data.items()))

# now write the sdt_data into a text file
filename_sdt= 'distances_' + folder_name + '.txt'
filepath_g_sdt = os.path.join(folder_name, filename_sdt)
with open(filepath_g_sdt, 'w') as file:
    for key, value in sdt_data.items():
        file.write(f'{key} {value}\n')







# # save boundary_mesh_polygon_n into a text file according to map_for_indices
# with open('boundary_mesh_polygon_n.txt', 'w') as file:
#     for vertex in boundary_mesh_polygon_n:
#         index1=vertices.index(vertex)
#         index2=map_for_indices[index1]
#         file.write(f'{index2}\n')

# #save boundary_mesh_polygon_n for every mesh_polygon_n into a text file according to map_for_indices 
# for i, poly in enumerate(mesh_polygons):
#     with open(f'boundary_mesh_polygon_{i}.txt', 'w') as file:
#         for vertex in eval(f'boundary_mesh_polygon_{i}'):
#             index1=vertices.index(vertex)
#             index2=eval(f'map_for_indices_{i}')[index1]
#             file.write(f'{index2}\n')


# sdt_file_path = 'blob0.sdt'
# sdt_data= read_sdt_file(sdt_file_path)


# # boundary_mesh_polygon_n_sdt.txt will contain the sdt data of boundary_mesh_polygon_n
# with open('boundary_mesh_polygon_n_sdt.txt', 'w') as file:
#     for vertex in boundary_mesh_polygon_n:
#         index1=vertices.index(vertex)
#         file.write(f'{sdt_data[index1]}\n')

# # boundary_mesh_polygon_n_sdt.txt will contain the sdt data of boundary_mesh_polygon_n for every mesh_polygon_n
# for i, poly in enumerate(mesh_polygons):
#     with open(f'boundary_mesh_polygon_{i}_sdt.txt', 'w') as file:
#         for vertex in eval(f'boundary_mesh_polygon_{i}'):
#             index1=vertices.index(vertex)
#             file.write(f'{sdt_data[index1]}\n')


# #plot the boundary mesh polygon n as scatter plot
# fig, ax = plt.subplots()
# x, y = zip(*boundary_mesh_polygon_n)
# ax.scatter(x, y, color='green', label='Boundary Points')
# plt.legend()
# plt.show()


# # Plot the polygons and mesh polygons
# plot_polygons(polygons, mesh_polygons)

#write an obj file for mesh_polygon_n using map_for_indices dictionary
# with open('mesh_polygon_n.obj', 'w') as file:
#     for i, vertex in enumerate(mesh_polygons[n]):
#         file.write(f"v {vertex[0]} {vertex[1]} 0\n")
#     for face in mesh_polygons_faces:
#         file.write(f"f {map_for_indices[face[0]]} {map_for_indices[face[1]]} {map_for_indices[face[2]]}\n")

# # write an obj file for every mesh_polygon_n using map_for_indices dictionary
# for i, poly in enumerate(mesh_polygons):
#     with open(f'mesh_polygon_{i}.obj', 'w') as file:
#         for j, vertex in enumerate(poly):
#             file.write(f"v {vertex[0]} {vertex[1]} 0\n")
#         for face in eval(f'mesh_polygons_faces_{i}'):
#             file.write(f"f {eval(f'map_for_indices_{i}')[face[0]]} {eval(f'map_for_indices_{i}')[face[1]]} {eval(f'map_for_indices_{i}')[face[2]]}\n")






















