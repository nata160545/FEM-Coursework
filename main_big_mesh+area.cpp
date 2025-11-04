#include "mesh.hpp"
#include <iostream>
#include <vector>

int main()
{
    using Value = float;
    using Mesh = fem::Mesh<Value>;
    using Size = Mesh::Size;
    // Большая сетка на всей области [0,1]x[0,1]
    typename Mesh::Nodes bigMeshNodes = Mesh::buildRegulArea(0.0, 0.0, 1.0, 1.0, 2, 2);
    Mesh bigMesh(bigMeshNodes);

    // Характестика большой сетки
    std::cout << "Big mesh created successfully!" << std::endl;
    std::cout << "Number of nodes: " << bigMeshNodes.size() << std::endl;
    std::cout << "Number of elements: " << bigMesh.getNumElements() << std::endl;

    // координаты узлов большой сетки
    std::cout << "\nNode coordinates of big mesh:" << std::endl;
    for (Size i = 0; i < bigMeshNodes.size(); ++i)
    {
        std::cout << "Node " << i << ": ("
                  << bigMeshNodes(i).coords(0) << ", "
                  << bigMeshNodes(i).coords(1) << ")" << std::endl;
    }

    // сохраняем в vtk
    std::cout << "\nSaving big mesh to VTK file..." << std::endl;
    bigMesh.writeParaViewVtk("big_mesh.vtk");
    std::cout << "Big mesh saved to 'big_mesh.vtk'" << std::endl;

    // 4 маленькие области
    typename Mesh::Nodes area0 = Mesh::buildRegulArea(0.0, 0.0, 0.5, 0.5, 2, 2);
    typename Mesh::Nodes area1 = Mesh::buildRegulArea(0.5, 0.0, 1.0, 0.5, 2, 2);
    typename Mesh::Nodes area2 = Mesh::buildRegulArea(0.0, 0.5, 0.5, 1.0, 2, 2);
    typename Mesh::Nodes area3 = Mesh::buildRegulArea(0.5, 0.5, 1.0, 1.0, 2, 2);

    Mesh mesh0(area0);
    Mesh mesh1(area1);
    Mesh mesh2(area2);
    Mesh mesh3(area3);

    // Сохраняем маленькие области
    mesh0.writeParaViewVtk("area0.vtk");
    mesh1.writeParaViewVtk("area1.vtk");
    mesh2.writeParaViewVtk("area2.vtk");
    mesh3.writeParaViewVtk("area3.vtk");

    std::cout << "- area0.vtk" << std::endl;
    std::cout << "- area1.vtk" << std::endl;
    std::cout << "- area2.vtk" << std::endl;
    std::cout << "- area3.vtk" << std::endl;
    return 0;
}
