#include "mesh.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

int main()
{
    using Value = float;
    using Mesh = fem::Mesh<Value>;
    using Size = Mesh::Size;

    // Большая сетка
    typename Mesh::Nodes bigMeshNodes = Mesh::buildRegulArea(0.0, 0.0, 1.0, 1.0, 2, 2);

    // Граничные условия
    for (Size i = 0; i < 3; ++i)
    { // нижние узлы
        bigMeshNodes(i).disps.resize(2);
        bigMeshNodes(i).disps(0) = {0, 0.0};
        bigMeshNodes(i).disps(1) = {1, 0.0};
    }
    for (Size i = 6; i < 9; ++i)
    { // верхние узлы
        bigMeshNodes(i).disps.resize(1);
        bigMeshNodes(i).disps(0) = {1, -0.1};
    }

    Mesh bigMesh(bigMeshNodes);
    Value elastMod = 200000;
    Value poissRat = 0.3;
    bigMesh.calculateStiffnessMatrix(elastMod, poissRat);
    bigMesh.calculateForceVector();
    bigMesh.calculateDisplacementVector();
    auto displacements = bigMesh.getDisplacementVector();
    std::cout << "\nBIG MESH DISPLACEMENTS:" << std::endl;
    std::cout << "Node\tu\t\tv" << std::endl;
    std::cout << "----------------------------" << std::endl;
    for (Size i = 0; i < bigMesh.getNumNodes(); ++i)
    {
        std::cout << i << "\t" << displacements(2 * i)
                  << "\t\t" << displacements(2 * i + 1) << std::endl;
    }
    bigMesh.writeParaViewVtk("big_mesh.vtk");

    // 4 маленькие области
    typename Mesh::Nodes area0 = Mesh::buildRegulArea(0.0, 0.0, 0.5, 0.5, 2, 2);
    typename Mesh::Nodes area1 = Mesh::buildRegulArea(0.5, 0.0, 1.0, 0.5, 2, 2);
    typename Mesh::Nodes area2 = Mesh::buildRegulArea(0.0, 0.5, 0.5, 1.0, 2, 2);
    typename Mesh::Nodes area3 = Mesh::buildRegulArea(0.5, 0.5, 1.0, 1.0, 2, 2);

    // выбираем область
    auto &targetArea = area1;

    for (Size i = 0; i < targetArea.size(); ++i)
    {
        targetArea(i).disps.resize(2);
        targetArea(i).disps(0) = {0, displacements(2 * i)};     // u
        targetArea(i).disps(1) = {1, displacements(2 * i + 1)}; // v

        std::cout << "node tsrgetarea " << i << ": u="
                  << displacements(2 * i) << ", v=" << displacements(2 * i + 1) << std::endl;
    }

    Mesh mesh0(area0);
    Mesh mesh1(area1);
    Mesh mesh2(area2);
    Mesh mesh3(area3);

    Mesh *targetMesh = nullptr;

    if (&targetArea == &area0)
        targetMesh = &mesh0;
    else if (&targetArea == &area1)
        targetMesh = &mesh1;
    else if (&targetArea == &area2)
        targetMesh = &mesh2;
    else if (&targetArea == &area3)
        targetMesh = &mesh3;

    mesh0.calculateStiffnessMatrix(elastMod, poissRat);
    mesh1.calculateStiffnessMatrix(elastMod, poissRat);
    mesh2.calculateStiffnessMatrix(elastMod, poissRat);
    mesh3.calculateStiffnessMatrix(elastMod, poissRat);

    // ВЫВОД МАТРИЦ ЖЁСТКОСТИ В ФАЙЛ
    std::ofstream stiffnessFile("stiffness_matrix.txt");
    if (stiffnessFile.is_open())
    {
        // Большая сетка
        stiffnessFile << "BIG MESH STIFFNESS MATRIX ("
                      << bigMesh.getStiffnessMatrix().rows() << "x"
                      << bigMesh.getStiffnessMatrix().cols() << "):\n";
        stiffnessFile << bigMesh.getStiffnessMatrix() << "\n\n";

        // Area 0
        stiffnessFile << "AREA0 STIFFNESS MATRIX ("
                      << mesh0.getStiffnessMatrix().rows() << "x"
                      << mesh0.getStiffnessMatrix().cols() << "):\n";
        stiffnessFile << mesh0.getStiffnessMatrix() << "\n\n";

        // Area 1
        stiffnessFile << "AREA1 STIFFNESS MATRIX ("
                      << mesh1.getStiffnessMatrix().rows() << "x"
                      << mesh1.getStiffnessMatrix().cols() << "):\n";
        stiffnessFile << mesh1.getStiffnessMatrix() << "\n\n";

        // Area 2
        stiffnessFile << "AREA2 STIFFNESS MATRIX ("
                      << mesh2.getStiffnessMatrix().rows() << "x"
                      << mesh2.getStiffnessMatrix().cols() << "):\n";
        stiffnessFile << mesh2.getStiffnessMatrix() << "\n\n";

        // Area 3
        stiffnessFile << "AREA3 STIFFNESS MATRIX ("
                      << mesh3.getStiffnessMatrix().rows() << "x"
                      << mesh3.getStiffnessMatrix().cols() << "):\n";
        stiffnessFile << mesh3.getStiffnessMatrix() << "\n";

        stiffnessFile.close();
        std::cout << "Stiffness matrix saved to 'stiffness_matrix.txt'" << std::endl;
    }
    else
    {
        std::cerr << "Error: Could not open stiffness_matrix.txt for writing!" << std::endl;
    }

    targetMesh->calculateStiffnessMatrix(elastMod, poissRat);
    targetMesh->calculateForceVector();
    targetMesh->calculateDisplacementVector();
    auto target_displacements = targetMesh->getDisplacementVector();
    std::cout << "\ndisplacement targetarea:" << std::endl;
    std::cout << "Node\tu\t\tv" << std::endl;
    std::cout << "----------------------------" << std::endl;
    for (Size i = 0; i < targetMesh->getNumNodes(); ++i)
    {
        std::cout << i << "\t" << target_displacements(2 * i)
                  << "\t\t" << target_displacements(2 * i + 1) << std::endl;
    }

    // ВЫВОД ПЕРЕМЕЩЕНИЙ В ТЕКСТОВЫЙ ФАЙЛ
    std::ofstream dispFile("displacements.txt");
    if (dispFile.is_open())
    {
        dispFile << std::fixed << std::setprecision(6);

        // Перемещения большой сетки
        dispFile << "BIG MESH DISPLACEMENTS:\n";
        dispFile << "Node\tX-Coord\t\tY-Coord\t\tu\t\tv\n";
        dispFile << "------------------------------------------------------------\n";
        for (Size i = 0; i < bigMesh.getNumNodes(); ++i)
        {
            dispFile << i << "\t"
                     << bigMeshNodes(i).coords(0) << "\t\t"
                     << bigMeshNodes(i).coords(1) << "\t\t"
                     << displacements(2 * i) << "\t\t"
                     << displacements(2 * i + 1) << "\n";
        }

        dispFile << "\n\nTARGET AREA DISPLACEMENTS:\n";
        dispFile << "Node\tX-Coord\t\tY-Coord\t\tu\t\tv\n";
        dispFile << "------------------------------------------------------------\n";
        for (Size i = 0; i < targetMesh->getNumNodes(); ++i)
        {
            dispFile << i << "\t"
                     << targetArea(i).coords(0) << "\t\t"
                     << targetArea(i).coords(1) << "\t\t"
                     << target_displacements(2 * i) << "\t\t"
                     << target_displacements(2 * i + 1) << "\n";
        }
    }
    mesh0.writeParaViewVtk("area0.vtk");
    mesh1.writeParaViewVtk("area1.vtk");
    mesh2.writeParaViewVtk("area2.vtk");
    mesh3.writeParaViewVtk("area3.vtk");

    if (&targetArea == &area0)
        targetMesh->writeParaViewVtk("target_area0.vtk");
    else if (&targetArea == &area1)
        targetMesh->writeParaViewVtk("target_area1.vtk");
    else if (&targetArea == &area2)
        targetMesh->writeParaViewVtk("target_area2.vtk");
    else if (&targetArea == &area3)
        targetMesh->writeParaViewVtk("target_area3.vtk");

    std::cout << "- big_mesh.vtk" << std::endl;
    std::cout << "- area0.vtk" << std::endl;
    std::cout << "- area1.vtk" << std::endl;
    std::cout << "- area2.vtk" << std::endl;
    std::cout << "- area3.vtk" << std::endl;

    if (&targetArea == &area0)
        std::cout << "- target_area0.vtk (displacements)" << std::endl;
    else if (&targetArea == &area1)
        std::cout << "- target_area1.vtk (displacements)" << std::endl;
    else if (&targetArea == &area2)
        std::cout << "- target_area2.vtk (displacements)" << std::endl;
    else if (&targetArea == &area3)
        std::cout << "- target_area3.vtk (displacements)" << std::endl;

    std::system("pause");
    return 0;
}
