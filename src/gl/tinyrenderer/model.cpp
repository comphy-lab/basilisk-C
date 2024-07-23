#include <iostream>
#include <sstream>
#include "geometry.h"
#include "model.h"

Model::Model(const std::string filename) {
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            vec3 v;
            iss >> v.x, iss >> v.y, iss >> v.z;
            verts.push_back(v);
        } else if (!line.compare(0, 3, "vn ")) {
            iss >> trash >> trash;
            vec3 n;
	    iss >> n.x, iss >> n.y, iss >> n.z;
            norms.push_back(vec3_normalized(n));
        } else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            vec2 uv;
	    iss >> uv.x, iss >> uv.y;
            tex_coord.push_back({uv.x, 1-uv.y});
        }  else if (!line.compare(0, 2, "f ")) {
            int f,t,n;
            iss >> trash;
            int cnt = 0;
            while (iss >> f >> trash >> t >> trash >> n) {
                facet_vrt.push_back(--f);
                facet_tex.push_back(--t);
                facet_nrm.push_back(--n);
                cnt++;
            }
            if (3!=cnt) {
                std::cerr << "Error: the obj file is supposed to be triangulated" << std::endl;
                return;
            }
        }
    }
    std::cerr << "# v# " << nverts() << " f# "  << nfaces() << " vt# " << tex_coord.size() << " vn# " << norms.size() << std::endl;
    load_texture(filename, "_diffuse.tga",    diffusemap );
    load_texture(filename, "_nm_tangent.tga", normalmap  );
    load_texture(filename, "_spec.tga",       specularmap);
}

int Model::nverts() const {
    return verts.size();
}

int Model::nfaces() const {
    return facet_vrt.size()/3;
}

vec3 Model::vert(const int i) const {
    return verts[i];
}

vec3 Model::vert(const int iface, const int nthvert) const {
    return verts[facet_vrt[iface*3+nthvert]];
}

void Model::load_texture(std::string filename, const std::string suffix, TGAImage &img) {
    size_t dot = filename.find_last_of(".");
    if (dot==std::string::npos) return;
    std::string texfile = filename.substr(0,dot) + suffix;
    std::cerr << "texture file " << texfile << " loading " << (img.read_tga_file(texfile.c_str()) ? "ok" : "failed") << std::endl;
}

vec3 Model::normal(const vec2 &uvf) const {
    TGAColor c = normalmap.get(uvf.x*normalmap.width(), uvf.y*normalmap.height());
    return vec3_sub ({c[2]*2./255.,c[1]*2./255.,c[0]*2./255.}, {1,1,1});
}

vec2 Model::uv(const int iface, const int nthvert) const {
    return tex_coord[facet_tex[iface*3+nthvert]];
}

vec3 Model::normal(const int iface, const int nthvert) const {
    return norms[facet_nrm[iface*3+nthvert]];
}
