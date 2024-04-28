def parse_obj_file(obj_file_path):
    vertices = []
    faces = []

    with open(obj_file_path, 'r') as obj_file:
        for line in obj_file:
            if line.startswith('v '):
                # Parse vertex coordinates
                _, x, y, z = line.strip().split()
                vertices.append((float(x), float(y), float(z)))
            elif line.startswith('f '):
                # Parse face vertices indices
                vertices_indices = [int(index.split('/')[0]) for index in line.strip().split()[1:]]
                # Subtract 1 from each index to convert to 0-based indexing
                faces.append(vertices_indices)

    return vertices, faces

def write_xml_output(vertices, faces, output_file_path):
    with open(output_file_path, 'w') as output_file:
        output_file.write('<vertexdata>\n')
        for vertex in vertices:
            output_file.write(f'\t{vertex[0]} {vertex[1]} {vertex[2]}\n')
        output_file.write('</vertexdata>\n\n')
        output_file.write('<faces>\n')
        for face in faces:
            if len(face) == 3:
                output_file.write(f'\t{" ".join(str(index - 1) for index in face)}\n')  # Convert to 0-based indexing
            elif len(face) > 3:
                # Split faces with more than 3 vertices into triangles
                for i in range(1, len(face) - 1):
                    output_file.write(f'\t{face[0] - 1} {face[i] - 1} {face[i + 1] - 1}\n')  # Convert to 0-based indexing
        output_file.write('</faces>\n')

def main():
    obj_file_path = 'Cat_v1_l3.obj'
    output_file_path = 'output.xml'

    vertices, faces = parse_obj_file(obj_file_path)
    write_xml_output(vertices, faces, output_file_path)

if __name__ == "__main__":
    main()
