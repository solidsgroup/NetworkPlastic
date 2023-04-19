/* MyGAL
 * Copyright (C) 2019 Pierre Vigier
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

 // STL
#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <fstream>
#include <math.h>
#include <numeric>  
#include <algorithm> 
// SFML
#include <SFML/Graphics.hpp>
// MyGAL
#include "FortuneAlgorithm.h"
#include "Vector2.h"

using namespace mygal;

using Float = double;

constexpr Float WindowWidth = 600.0f*2.0f;
constexpr Float WindowHeight = 600.0f*2.0f;
constexpr Float PointRadius = 0.005f;
constexpr Float Offset = 1.0f;

// Points generation

template<typename T>
std::vector<Vector2<T>> generatePoints(int nbPoints)
{
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::cout << "seed: " << seed << '\n';
    auto generator = std::default_random_engine(seed);
    auto distribution = std::uniform_real_distribution<T>(0.0, 1.0);

    auto points = std::vector<Vector2<T>>(nbPoints);
    for (auto i = 0; i < nbPoints; ++i)
        points[i] = Vector2<T>(distribution(generator), distribution(generator));

    return points;
}
template<typename T>
std::vector<Vector2<T>> generateBimodalLaminate(int numSec, int finePnts, int coarsePnts)
{
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::cout << "seed: " << seed << '\n';
    auto generator = std::default_random_engine(seed);
    auto distribution_x = std::uniform_real_distribution<T>(0.0, 1.0);
    std::vector<Vector2<T>> points;
    for (int i = 0; i < numSec; ++i) {
        double final = 1.0/numSec + double(1.0 / numSec*i);
        double start = 0.0 + double(1.0/numSec*i);
        // std::cout << "start = " << start << " " << " finish = " << final << std::endl;
        auto distribution_y = std::uniform_real_distribution<T>(start, final);
        if (i % 2 == 1) {
            for (int j = 0; j < coarsePnts; j++) {
                points.push_back(Vector2<T>(distribution_x(generator), distribution_y(generator)));
            }
        }
        else {
            for (int j = 0; j < finePnts; j++) {
                points.push_back(Vector2<T>(distribution_x(generator), distribution_y(generator)));
            }
        }
    }

    return points;
}

// Rendering

template<typename T>
void drawPoint(sf::RenderWindow& window, Vector2<T> point, sf::Color color)
{
    auto shape = sf::CircleShape(PointRadius);
    shape.setPosition(sf::Vector2f(point.x - PointRadius, 1.0 - point.y - PointRadius));
    shape.setFillColor(color);
    window.draw(shape);
}

template<typename T>
void drawEdge(sf::RenderWindow& window, Vector2<T> origin, Vector2<T> destination, sf::Color color)
{
    auto line = std::array<sf::Vertex, 2>
    {
        sf::Vertex(sf::Vector2f(origin.x, 1.0 - origin.y), color),
            sf::Vertex(sf::Vector2f(destination.x, 1.0 - destination.y), color)
    };
    window.draw(line.data(), 2, sf::Lines);
}

template<typename T>
void drawPoints(sf::RenderWindow& window, const Diagram<T>& diagram)
{
    for (const auto& site : diagram.getSites())
        drawPoint(window, site.point, sf::Color(100, 250, 50));
}

template<typename T>
void drawDiagram(sf::RenderWindow& window, const Diagram<T>& diagram)
{
    for (const auto& site : diagram.getSites())
    {
        auto center = site.point;
        auto face = site.face;
        auto halfEdge = face->outerComponent;
        if (halfEdge == nullptr)
            continue;
        while (halfEdge->prev != nullptr)
        {
            halfEdge = halfEdge->prev;
            if (halfEdge == face->outerComponent)
                break;
        }
        auto start = halfEdge;
        while (halfEdge != nullptr)
        {
            if (halfEdge->origin != nullptr && halfEdge->destination != nullptr)
            {
                auto origin = (halfEdge->origin->point - center) * Offset + center;
                auto destination = (halfEdge->destination->point - center) * Offset + center;
                drawEdge(window, origin, destination, sf::Color::Blue);
            }
            halfEdge = halfEdge->next;
            if (halfEdge == start)
                break;
        }
    }
}

template<typename T>
void drawTriangulation(sf::RenderWindow& window, const Diagram<T>& diagram, const Triangulation& triangulation)
{
    for (auto i = std::size_t(0); i < diagram.getNbSites(); ++i)
    {
        auto origin = diagram.getSite(i)->point;
        for (const auto& j : triangulation.getNeighbors(i))
        {
            auto destination = diagram.getSite(j)->point;
            drawEdge(window, origin, destination, sf::Color::Green);
        }
    }
}

// Generating the diagram

template<typename T>
Diagram<T> generateDiagram(const std::vector<Vector2<T>>& points)
{
    // Construct diagram
    auto algorithm = FortuneAlgorithm<T>(points);
    auto start = std::chrono::steady_clock::now();
    algorithm.construct();
    auto duration = std::chrono::steady_clock::now() - start;
    std::cout << "construction: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "ms" << '\n';

    // Bound the diagram
    start = std::chrono::steady_clock::now();
    algorithm.bound(Box<T>{-0.05, -0.05, 1.05, 1.05}); // Take the bounding box slightly bigger than the intersection box
    duration = std::chrono::steady_clock::now() - start;
    std::cout << "bounding: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "ms" << '\n';
    auto diagram = algorithm.getDiagram();

    // Intersect the diagram with a box
    start = std::chrono::steady_clock::now();
    diagram.intersect(Box<T>{0.0, 0.0, 1.0, 1.0});
    duration = std::chrono::steady_clock::now() - start;
    std::cout << "intersection: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << "ms" << '\n';

    return diagram;
}

struct graph
{
    std::vector<int> s;
    std::vector<int> t;
    int edges;
    int N;
    std::vector<std::vector<mygal::Vector2<double>>> vertices;
    std::vector<double> area;
    std::vector<double> angle;
    std::vector<double> GBarea;
    std::vector<bool> flagGrain;
};

template<typename T>
graph MakeGraph(const Diagram<T>& dia,int pnts, std::vector<mygal::Vector2<double>> points,bool flag)
{
    graph G;
    std::vector<int> tempS;
    std::vector<int> tempT;
    std::vector<bool> del;
    G.area = dia.computeArea(); // get area of each polygon
    std::vector<std::string> flagGrain;
    G.vertices = dia.ReturnV(); // get vertices of each polygon
    auto triangulation = dia.computeTriangulation();  // construct graph
    double tol = 1e-20;

    for (int i = 0; i < pnts; i++) {
        auto& Neighbor = triangulation.getNeighbors(i);
        G.flagGrain.emplace_back(flag);
        for (int j = 0; j < Neighbor.size(); j++) {
            tempS.push_back(i);
            tempT.push_back(Neighbor.at(j));
            del.push_back(false);
        }
    }
    // delete duplicates, each edge has duplicate inverse value
    for (int i = 0; i < tempS.size(); i++) {
        int vert_s = tempS.at(i);
        int vert_t = tempT.at(i);
        for (int j = 0; j < tempT.size(); j++) {
            if (i == j ) continue;
            if (del.at(i)) break;
            int test_s = tempT.at(j);
            int test_t = tempS.at(j);
            if (vert_s == test_s &&
                vert_t == test_t) {
                del.at(j) = true;
                break;
            }
        }
    }
    for (int i = 0; i < tempS.size(); i++) {
        if (!del.at(i)) {
            G.s.push_back(tempS.at(i));
            G.t.push_back(tempT.at(i));
        }
    }

    // compute a_pq and inclination of boundary
    std::vector<mygal::Vector2<double>> triPnt;
    for (int i = 0; i < G.s.size(); i++) {
        std::vector<mygal::Vector2<double>> Vs = G.vertices.at(G.s.at(i));
        std::vector<mygal::Vector2<double>> Vt = G.vertices.at(G.t.at(i));
        triPnt.clear();
        for (int j = 0; j < Vs.size(); j++) {
            for (int k = 0; k < Vt.size(); k++) {
                if ((Vs.at(j)-Vt.at(k)).getNorm() <= tol) {
                    triPnt.push_back(Vs.at(j));
                    if (triPnt.size() == 2) {
                        G.GBarea.push_back((triPnt.at(0) - triPnt.at(1)).getNorm());
                        G.angle.push_back(triPnt.at(1).getAngle(triPnt.at(0)));
                        break;
                    }
                }
            }
        }
    }
    G.edges = G.s.size();
    G.N = G.edges + pnts;
    return G;
}


struct loc
{
    std::vector<int> grain_id_top;
    std::vector<int> grain_id_bottom;
    std::vector<double> gbArea_top;
    std::vector<double> gbArea_bottom;
    std::vector<double> start_vert_top;
    std::vector<double> end_vert_top;
    std::vector<double> start_vert_bottom;
    std::vector<double> end_vert_bottom;
};
graph MergeG(std::vector<graph>& G) {
    double tol = 1e-15; 

    std::vector<loc> id;
    graph MG;
    MG.edges = 0;
    MG.N = 0;
    for (int Gnum = 0; Gnum < G.size(); Gnum++) { // get top or bottom grains
        loc idtemp;
        std::vector<std::vector<mygal::Vector2<double>>> verts = G.at(Gnum).vertices;
        for (int vertex = 0; vertex < G.at(Gnum).area.size(); vertex++) {
            std::vector<mygal::Vector2<double>> pnts = verts.at(vertex);
            std::vector<mygal::Vector2<double>> pnts_temp1;
            std::vector<mygal::Vector2<double>> pnts_temp2;
            int grain_id;
            bool add_top = false;
            bool add_bottom = false;
            for (int i = 0; i < pnts.size(); i++) {
                if (std::abs(pnts.at(i).y - 1.0) <= tol) { // id top
                    pnts_temp1.emplace_back(pnts.at(i));
                    grain_id = vertex;
                    add_top = true;
                }
                else if (std::abs(pnts.at(i).y) <= tol) { // id bottom
                    pnts_temp2.emplace_back(pnts.at(i));
                    grain_id = vertex;
                    add_bottom = true;
                }
                G.at(Gnum).vertices.at(vertex).at(i).y += double(Gnum);
            }
            if (add_top) {
                idtemp.grain_id_top.emplace_back(grain_id);
                idtemp.gbArea_top.emplace_back(std::abs(pnts_temp1.at(0).x - pnts_temp1.at(1).x));
                idtemp.start_vert_top.emplace_back(std::min(pnts_temp1.at(0).x, pnts_temp1.at(1).x));
                idtemp.end_vert_top.emplace_back(std::max(pnts_temp1.at(0).x, pnts_temp1.at(1).x));
            }
            if (add_bottom) {
                idtemp.grain_id_bottom.emplace_back(grain_id);
                idtemp.gbArea_bottom.emplace_back(std::abs(pnts_temp2.at(0).x - pnts_temp2.at(1).x));
                idtemp.start_vert_bottom.emplace_back(std::min(pnts_temp2.at(0).x, pnts_temp2.at(1).x));
                idtemp.end_vert_bottom.emplace_back(std::max(pnts_temp2.at(0).x, pnts_temp2.at(1).x));
            }
        }
        id.emplace_back(idtemp);
    }
    for (int Gnum = 0; Gnum < G.size(); Gnum++) { // order grains from x=0 -> x=1
        double start = 0.0;
        loc idtemp; 
        loc idhold = id.at(Gnum);
        int i = 0; 
        while(start - 1.0 <= tol) {
            if (std::abs(idhold.start_vert_top.at(i) - start) <= tol) {
                idtemp.grain_id_top.emplace_back(idhold.grain_id_top.at(i));
                idtemp.gbArea_top.emplace_back(idhold.gbArea_top.at(i));
                idtemp.start_vert_top.emplace_back(idhold.start_vert_top.at(i));
                idtemp.end_vert_top.emplace_back(idhold.end_vert_top.at(i));

                start = idhold.end_vert_top.at(i);

                idhold.end_vert_top.erase(idhold.end_vert_top.begin() + i);
                idhold.grain_id_top.erase(idhold.grain_id_top.begin() + i);
                idhold.gbArea_top.erase(idhold.gbArea_top.begin() + i);
                idhold.start_vert_top.erase(idhold.start_vert_top.begin() + i);
            }
            i++;
            if (idhold.start_vert_top.size() == 0) {
                break;
            }
            i = i % idhold.start_vert_top.size();
        }
        start = 0.0;
        i = 0;
        while (start - 1.0 <= tol) {
            if (std::abs(idhold.start_vert_bottom.at(i) - start) <= tol) {
                idtemp.grain_id_bottom.emplace_back(idhold.grain_id_bottom.at(i));
                idtemp.gbArea_bottom.emplace_back(idhold.gbArea_bottom.at(i));
                idtemp.start_vert_bottom.emplace_back(idhold.start_vert_bottom.at(i));
                idtemp.end_vert_bottom.emplace_back(idhold.end_vert_bottom.at(i));

                start = idhold.end_vert_bottom.at(i);

                idhold.end_vert_bottom.erase(idhold.end_vert_bottom.begin() + i);
                idhold.grain_id_bottom.erase(idhold.grain_id_bottom.begin() + i);
                idhold.gbArea_bottom.erase(idhold.gbArea_bottom.begin() + i);
                idhold.start_vert_bottom.erase(idhold.start_vert_bottom.begin() + i);
            }
            i++;
            if (idhold.start_vert_bottom.size() == 0) {
                id.at(Gnum) = idtemp;
                break;
            }
            i = i % idhold.start_vert_bottom.size();
        }


    }
    int t_offset = 0;
    int s_offset = 0;
    for (int Gnum = 0; Gnum < G.size(); Gnum++) { // add values from individual graphs to one 
        
        MG.N += G.at(Gnum).area.size();
        s_offset = t_offset;
        for (int vertex = 0; vertex < G.at(Gnum).area.size(); vertex++) {
            
            MG.area.emplace_back(G.at(Gnum).area.at(vertex));
            MG.vertices.emplace_back(G.at(Gnum).vertices.at(vertex));
            MG.flagGrain.emplace_back(G.at(Gnum).flagGrain.at(vertex));
        }
        for (int edge = 0; edge < G.at(Gnum).s.size(); edge++) {

            MG.s.emplace_back(G.at(Gnum).s.at(edge) + s_offset);
            MG.t.emplace_back(G.at(Gnum).t.at(edge) + t_offset);
            MG.GBarea.emplace_back(G.at(Gnum).GBarea.at(edge));
            MG.angle.emplace_back(G.at(Gnum).angle.at(edge));
        }
        t_offset += G.at(Gnum).N - G.at(Gnum).edges;
    }
    t_offset = 0;
    s_offset = 0;
    for (int gnum = 0; gnum < id.size()-1; gnum++) { // connect top and bottom of diagrams
        t_offset += G.at(gnum).N - G.at(gnum).edges;
        std::vector<int> coarse_id;
        std::vector<double> coarse_gbArea;
        std::vector<int> fine_id;
        std::vector<double> fine_gbArea;
        if (gnum % 2 == 0) {
            coarse_id = id.at(gnum + 1).grain_id_bottom;
            coarse_gbArea = id.at(gnum + 1).gbArea_bottom;

            fine_id = id.at(gnum).grain_id_top;
            fine_gbArea = id.at(gnum).gbArea_top;
        }
        else {
            coarse_id = id.at(gnum).grain_id_top;
            coarse_gbArea = id.at(gnum).gbArea_top;

            fine_id = id.at(gnum + 1).grain_id_bottom;
            fine_gbArea = id.at(gnum + 1).gbArea_bottom;
        }
        int j = 0; 
        double bottom_diff = 0.0;
        double top_diff = 0.0;
        for (int i = 0; i < coarse_gbArea.size(); i++) { // coarse grains

            bottom_diff += coarse_gbArea.at(i);
            double total = 0.0;
            for (j; j < fine_gbArea.size(); j++)  { // fine grains

                top_diff += fine_gbArea.at(j);
                total += fine_gbArea.at(j);
                if (top_diff <= bottom_diff) {
                    int s = fine_id.at(j) + s_offset;
                    int t = coarse_id.at(i) + t_offset;
                    MG.s.emplace_back(s);
                    MG.t.emplace_back(t);
                    MG.GBarea.emplace_back(fine_gbArea.at(j));
                    MG.angle.emplace_back(0.0);
                }
                else if (top_diff > bottom_diff) {
                    int s = fine_id.at(j) + s_offset;
                    int t = coarse_id.at(i) + t_offset;
                    MG.s.emplace_back(s);
                    MG.t.emplace_back(t);
                    MG.GBarea.emplace_back(coarse_gbArea.at(i) - (total - fine_gbArea.at(j)));
                    MG.angle.emplace_back(0.0);

                    s = fine_id.at(j) + s_offset;
                    t = coarse_id.at(i+1) + t_offset;
                    MG.s.emplace_back(s);
                    MG.t.emplace_back(t);
                    MG.GBarea.emplace_back(total - coarse_gbArea.at(i));
                    MG.angle.emplace_back(0.0);
                    j++;
                    break;
                }
            }
        }
        s_offset = s_offset;
    }
    MG.edges = MG.s.size();
    MG.N += MG.edges;
    return MG;
}
void saveMicrostructure(graph& G, bool bimodal,bool NP_Stat, int numSec, std::string fname = "Microstructure.txt")
{
    int NumVertices = G.N - G.edges;
    std::string path; std::string path_poly;
    if (NP_Stat) {
        path = "NP_Stat/" + fname;
        path_poly = "NP_Stat/poly_"  + fname;
    }

    else {
        path = fname;
        path_poly = "poly_" + fname;
    }

    std::ofstream myfile(path, std::ios::trunc);
    myfile << "s, t, a_{pq}, GB_Angle, Volume, Flag Grain, (number vertices = " << NumVertices << ")" << std::endl;
    myfile << NumVertices << " " << bimodal << " " << numSec << " " << NumVertices << " " << NumVertices << " " << NumVertices << std::endl;
    if (myfile.is_open())
    {
        for (int val = 0; val < G.edges; val++) {

            myfile << G.s.at(val) << " " << G.t.at(val) << " " << G.GBarea.at(val) << " " << G.angle.at(val); // edge values
            if (val < NumVertices) {

                myfile << " " << G.area.at(val) << " " << G.flagGrain.at(val); // vertex values
            }
            myfile << std::endl;
        }
        myfile.close();
    }
    else std::cout << "Unable to open microstructure file";

    // save vertices of grains for reconstruction
    std::ofstream myfile2(path_poly, std::ios::trunc);
    
    if (myfile2.is_open()) {

        for (int i = 0; i < NumVertices; i++) {
            auto& vert = G.vertices.at(i);
            for (int j = 0; j < vert.size(); j++) {

                myfile2 << vert.at(j).x << " " << vert.at(j).y << " ";

            }
            myfile2 << std::endl;
        }
        myfile2.close();
    }
}


int main()
{
    bool NP_Stat = true;
    int numSec = 8; int finePnts = 500; int coarsePnts = 3; int N = log2(numSec); // N1 - 11001,3 N2 - 9001, 4 N3 - 11001, 3

    bool bimodal = true;
    //int nbPoints = std::floor(numSec / 2) * coarsePnts + (numSec - std::floor(numSec / 2)) * finePnts; 
    //std::cout << "Number of grians = " << nbPoints << std::endl;
    
    for (int i = 0; i < 100; i++) { 
        //std::string fname = "microstructure_test_N"  + std::to_string(N) + "_" + std::to_string(i);
        std::string fname = "microstructure_small_" + std::to_string(i);
        std::vector<graph> Gs;
        std::cout << "iter = " << i << std::endl;

        auto bimodal_pnts = generatePoints<Float>(finePnts);
        auto diagram = generateDiagram(bimodal_pnts);
        diagram = generateDiagram(diagram.computeLloydRelaxation());
        auto triangulation = diagram.computeTriangulation();
        bool flag = false;
        graph G = MakeGraph(diagram, bimodal_pnts.size(), bimodal_pnts, flag);
        
        //for (int j = 0; j < numSec/2; j++) {
         //   auto bimodal_pnts = generatePoints<Float>(finePnts);
          //  auto diagram = generateDiagram(bimodal_pnts);
           // diagram = generateDiagram(diagram.computeLloydRelaxation());
           // auto triangulation = diagram.computeTriangulation();
           // bool flag = false;
           // graph G1 = MakeGraph(diagram, bimodal_pnts.size(), bimodal_pnts, flag);
            
          //  bimodal_pnts = generatePoints<Float>(coarsePnts);
          //  auto diagram2 = generateDiagram(bimodal_pnts);
          //  diagram2 = generateDiagram(diagram2.computeLloydRelaxation());
          //  auto triangulation2 = diagram2.computeTriangulation();
           // flag = true;
           // graph G2 = MakeGraph(diagram2, bimodal_pnts.size(), bimodal_pnts, flag);

          //  Gs.emplace_back(G1);
           // Gs.emplace_back(G2);
       // }
        //graph G = MergeG(Gs);
        saveMicrostructure(G, bimodal, NP_Stat, numSec, fname);   
    }
    
    auto bimodal_pnts = generatePoints<Float>(finePnts);
    auto diagram = generateDiagram(bimodal_pnts);
    diagram = generateDiagram(diagram.computeLloydRelaxation());
    auto triangulation = diagram.computeTriangulation();
    // Display the diagram
    auto settings = sf::ContextSettings();
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WindowWidth, WindowHeight), "MyGAL", sf::Style::Default, settings); // Can use auto only in C++17
    window.setView(sf::View(sf::FloatRect(-0.1f, -0.1f, 1.2f, 1.2f)));
    window.setFramerateLimit(60);
    auto showTriangulation = false;

    while (window.isOpen())
    {
        auto event = sf::Event();
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            else if (event.type == sf::Event::KeyReleased)
            {
                if (event.key.code == sf::Keyboard::Key::N)
                {
                    diagram = generateDiagram(generatePoints<Float>(bimodal_pnts.size()));
                    triangulation = diagram.computeTriangulation();
                }
                else if (event.key.code == sf::Keyboard::Key::R)
                {
                    diagram = generateDiagram(diagram.computeLloydRelaxation());
                    triangulation = diagram.computeTriangulation();
                }
                else if (event.key.code == sf::Keyboard::Key::T)
                    showTriangulation = !showTriangulation;
            }
        }
        window.clear(sf::Color::White);

        if (!showTriangulation) {
            drawDiagram(window, diagram);
            //drawPoints(window, diagram);
        }
        if (showTriangulation) {
            drawTriangulation(window, diagram, triangulation);
        }

        window.display();
    }

    
    return 0;
}
