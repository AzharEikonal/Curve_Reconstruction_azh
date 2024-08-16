import numpy as np
import matplotlib.pyplot as plt
import trimesh
from scipy.spatial import cKDTree
import sys
import os

def parse_obj_file(obj_file_path):
    vertices = []
    faces = []
    with open(obj_file_path, 'r') as file:
        for line in file:
            if line.startswith('v '):
                parts = line.split()
                vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('f '):
                parts = line.split()
                face = [int(part.split('/')[0]) for part in parts[1:]]
                faces.append(face)
    return vertices, faces

def transfinite_mean_value_interpolation(vertices, boundary_indices, boundary_distances):
    dist=[]
    for i in range(len(vertices)):
        for k in range(len(boundary_indices)):
            if i==boundary_indices[k]:
                dist.append(boundary_distances[k])
                break
        else:
            sum1=0
            sum2=0
            for j in range(len(boundary_indices)):
                d=np.linalg.norm(vertices[i]- vertices[boundary_indices[j]])
                quotient1= boundary_distances[j]/ d
                sum1+=quotient1
                quotient2= 1/d
                sum2+=quotient2
            dist.append(sum1/sum2)
    return dist 




# vertices, faces = parse_obj_file('mesh_polygon_n.obj')
# vertices = np.array(vertices)
# faces = np.array(faces)

if len(sys.argv) < 2:
    print("Please provide a folder name as an argument.")
    sys.exit(1)

folder_name = sys.argv[1]

# count the number of mesh_polygon_n files in the folder
count = 0
for file_name in os.listdir(folder_name):
    if file_name.startswith('mesh_polygon_') and file_name.endswith('.obj'):
        count += 1




# load vertices and faces for every mesh polygon_n in vertices_n and faces_n array 
for i in range(count):
    file_path_obj = os.path.join(folder_name, f'mesh_polygon_{i}.obj')  
    string = f'vertices_{i}, faces_{i} = parse_obj_file(file_path_obj)'
    exec(string)
    string = f'vertices_{i} = np.array(vertices_{i})'
    exec(string)
    string = f'faces_{i} = np.array(faces_{i})'
    exec(string)


# #load the bounary vertices from a txt file
# boundary_vertices = []
# with open('boundary_mesh_polygon_n.txt', 'r') as file:
#     for line in file:
#         boundary_vertices.append(int(line))
# boundary_vertices = np.array(boundary_vertices)

# load boundary vertices for every mesh polygon_n
for i in range(count):
    string =f'boundary_vertices_{i} = []'
    exec(string)
    file_path_txt = os.path.join(folder_name, f'boundary_mesh_polygon_{i}.txt')
    with open(file_path_txt, 'r') as file:
        for line in file:
            string = f'boundary_vertices_{i}.append(int(line))'
            exec(string)
    string = f'boundary_vertices_{i} = np.array(boundary_vertices_{i})'
    exec(string)




# boundary_distancs=[]
# with open('boundary_mesh_polygon_n_sdt.txt', 'r') as file:
#     for line in file:
#         boundary_distancs.append(float(line))
# boundary_distancs = np.array(boundary_distancs)

# load boundary distances for every mesh polygon_n
for i in range(count):
    string = f'boundary_distances_{i} = []'
    exec(string)
    file_path_dist_txt = os.path.join(folder_name, f'boundary_mesh_polygon_{i}_sdt.txt')
    with open(file_path_dist_txt, 'r') as file:
        for line in file:
            string = f'boundary_distances_{i}.append(float(line))'
            exec(string)
    string = f'boundary_distances_{i} = np.array(boundary_distances_{i})'
    exec(string)

# interpolate distances for every mesh polygon_n
for i in range(count):
    string= f'distances_{i} = transfinite_mean_value_interpolation(vertices_{i}, boundary_vertices_{i}, boundary_distances_{i})'
    exec(string)



# distances =transfinite_mean_value_interpolation(vertices, boundary_vertices, boundary_distancs)

# # save distances in a txt file
# with open('distances_mesh_polygon_n.txt', 'w') as file:
#     for distance in distances:
#         file.write(str(distance) + '\n')

# save distances for every mesh polygon_n in a txt file
for i in range(count):
    file_path_dist_txt = os.path.join(folder_name, f'distances_mesh_polygon_{i}.txt')
    with open(file_path_dist_txt, 'w') as file:
        string = f'for distance in distances_{i}:\n\tfile.write(str(distance) + \'\\n\')'
        exec(string)
    





