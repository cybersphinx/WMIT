/*
	Copyright 2010 Warzone 2100 Project

	This file is part of WMIT.

	WMIT is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	WMIT is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with WMIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Mesh.hpp"

#include <algorithm>
#include <iterator>
#include <map>
#include <set>

#include <cmath>

#include <sstream>

#include "Generic.hpp"
#include "Util.hpp"
#include "Pie.hpp"
#include "Vector.hpp"


WZMConnector::WZMConnector(GLfloat x, GLfloat y, GLfloat z):
	m_pos(Vertex<GLfloat>(x, y, z))
{
}

const WZMVertex& WZMConnector::getPos() const
{
	return m_pos;
}

Mesh::Mesh()
{
	defaultConstructor();
}

Mesh::Mesh(const Pie3Level& p3, float uvEps, float vertEps)
{
	std::vector<Pie3Polygon>::const_iterator itL;

	WZMVertex::less_wEps compare(vertEps);
	std::multimap<WZMVertex, unsigned, WZMVertex::less_wEps> tmpMMap(compare);

	std::multimap<WZMVertex, unsigned>::iterator itMM;
	std::pair<std::multimap<WZMVertex, unsigned>::iterator,
				std::multimap<WZMVertex, unsigned>::iterator> ret;

	unsigned vIdx;
	unsigned vert;
	unsigned frame;
	bool foundMatch;

	defaultConstructor();

	/*
	 *	Try to prevent duplicate vertices
	 *	(remember, different UV's, or animations,
	 *	 will cause unavoidable duplication)
	 *	so that our transformed vertex cache isn't
	 *	completely useless.
	 */

	if (m_textureArrays.size() == 0)
	{
		m_textureArrays.push_back(TexArray());
	}

	// For each pie3 polygon
	for (itL = p3.m_polygons.begin(); itL != p3.m_polygons.end(); ++itL)
	{
		IndexedTri iTri;

		// For all 3 vertices of the triangle
		for (vert = 0; vert < 3; ++vert)
		{
			WZMVertex wzmVert = p3.m_points[itL->getIndex(vert)];
			ret = tmpMMap.equal_range(wzmVert);
			foundMatch = false;

			// Only check if match is possible
			if (itL->m_frames <= m_textureArrays.size())
			{
				// For each potential match
				for (itMM = ret.first; itMM != ret.second; ++itMM)
				{
					vIdx = itMM->second;

					// Compare each animation frame
					for (frame = 0; frame < itL->m_frames; ++frame)
					{
						// Approximate comparison, helps kill duplicates
						const WZMUV::equal_wEps compare(uvEps);
						if (!compare(m_textureArrays[frame][vIdx], itL->getUV(vert, frame)))
						{
							break; // Not equal
						}
					}

					// if all frames were equal
					if (!(frame < itL->m_frames))
					{
						foundMatch = true;
						break; // Stop looking
					}
				}
			}

			if (!foundMatch)
			{
				unsigned frames2Do = std::max<size_t>(itL->m_frames, m_textureArrays.size());
				vIdx = m_vertexArray.size();
				// add the vertex to both the multimap and the vertex array
				m_vertexArray.push_back(wzmVert);
				tmpMMap.insert(std::make_pair(wzmVert, vIdx));

				m_textureArrays.reserve(frames2Do);

				// add the uv's to the texture arrays
				for (frame = 0; frame < frames2Do; ++frame)
				{
					if (m_textureArrays.size() < frame + 1)
					{
						// Expand the texture array arrays
						m_textureArrays.push_back(m_textureArrays.back());
						m_textureArrays[frame][vIdx] = itL->getUV(vert, frame % itL->m_frames);
					}
					else
					{
						m_textureArrays[frame].push_back(itL->getUV(vert, frame % itL->m_frames));
					}
				}
			}

			// Set the index
			iTri[vert] = vIdx;
		}

		m_indexArray.push_back(iTri);
	}

	std::list<Pie3Connector>::const_iterator itC;

	// For each pie3 connector
	for (itC = p3.m_connectors.begin(); itC != p3.m_connectors.end(); ++itC)
	{
		addConnector(WZMConnector(itC->pos.operator[](0),
                                  itC->pos.operator[](1),
                                  itC->pos.operator[](2)));
	}

	recalculateBoundData();
}

Mesh::~Mesh()
{
}

Pie3Level Mesh::backConvert(const Mesh& wzmMesh)
{
	return wzmMesh;
}

Mesh::operator Pie3Level() const
{
	Pie3Level p3;

	std::vector<Pie3Vertex>::iterator itPV;

	std::vector<IndexedTri>::const_iterator itTri;

	unsigned i;

	/* Note:
	 * WZM will have duplicates due to uv map differences
	 * so we remove those when converting
	 */

	for (itTri = m_indexArray.begin(); itTri != m_indexArray.end(); ++itTri)
	{
		Pie3Polygon p3Poly;
		Pie3UV	p3UV;

		p3Poly.m_flags = 0x200;

		for (i = 0; i < 3; ++i)
		{
			typedef Pie3Vertex::equal_wEps equals;
			mybinder1st<equals> compare(m_vertexArray[(*itTri)[i]]);

			itPV = std::find_if(p3.m_points.begin(), p3.m_points.end(), compare);

			if (itPV == p3.m_points.end())
			{
				// add it now
				p3Poly.m_indices[i] = p3.m_points.size();
				p3.m_points.push_back(m_vertexArray[(*itTri)[i]]);
			}
			else
			{
				p3Poly.m_indices[i] = std::distance(p3.m_points.begin(), itPV);
			}

			// TODO: deal with UV animation
			p3UV.u() = m_textureArrays[0][(*itTri)[i]].u();
			p3UV.v() = m_textureArrays[0][(*itTri)[i]].v();
			p3Poly.m_texCoords[i] = p3UV;
		}
		p3.m_polygons.push_back(p3Poly);
	}

	std::list<WZMConnector>::const_iterator itC;

	// For each WZM connector
	for (itC = m_connectors.begin(); itC != m_connectors.end(); ++itC)
	{
		Pie3Connector conn;
        conn.pos.operator[](0) = itC->getPos().operator[](0);
		conn.pos.operator[](1) = itC->getPos().operator[](1);
		conn.pos.operator[](2) = itC->getPos().operator[](2);
		p3.m_connectors.push_back(conn);
	}

	return p3;
}

bool Mesh::read(std::istream& in)
{
	std::string str;
	unsigned i,j,vertices,indices;
	GLfloat f;

	clear();

	in >> str >> m_name;
	if (in.fail() || str.compare("MESH") != 0)
	{
		std::cerr << "Mesh::read - Expected MESH directive found " << str;
		clear();
		return false;
	}

	if (!isValidWzName(m_name))
	{
		std::cerr << "Mesh::read - Invalid Mesh name.";
		m_name = std::string();
	}

	in >> str >> m_teamColours;
	if (in.fail() || str.compare("TEAMCOLOURS") != 0)
	{
		std::cerr << "Mesh::read - Expected TEAMCOLOURS directive found " << str;
		clear();
		return false;
	}

	in >> str >> vertices;
	if (in.fail() || str.compare("VERTICES") != 0)
	{
		std::cerr << "Mesh::read - Expected VERTICES directive found " << str;
		clear();
		return false;
	}

	in >> str >> indices;
	if (in.fail() || str.compare("FACES") != 0)
	{
		std::cerr << "Mesh::read - Expected FACES directive found " << str;
		clear();
		return false;
	}

	in >> str;
	if (str.compare("VERTEXARRAY") !=0)
	{
		std::cerr << "Mesh::read - Expected VERTEXARRAY directive found " << str;
		clear();
		return false;
	}

	m_vertexArray.reserve(vertices);
	for (; vertices > 0; --vertices)
	{
		WZMVertex vert;
		in >> vert.x() >> vert.y() >> vert.z();
		if (in.fail())
		{
			std::cerr << "Mesh::read - Error reading number of faces";
			clear();
			return false;
		}
		m_vertexArray.push_back(vert);
	}

	in >> str >> i;
	if (in.fail() || str.compare("TEXTUREARRAYS") != 0)
	{
		std::cerr << "Mesh::read - Expected TEXTUREARRAYS directive found " << str;
		clear();
		return false;
	}

	m_vertexArray.reserve(i);
	for (; i > 0; --i)
	{
		std::vector<WZMUV> tmp;
		tmp.clear();

		// j is currently ignored
		in >> str >> j;
		if ( in.fail() || str.compare("TEXTUREARRAY") != 0)
		{
			std::cerr << "Mesh::read - Expected TEXTUREARRAY directive found " << str;
			clear();
			return false;
		}

		for (j = 0; j < m_vertexArray.size(); ++j)
		{
			WZMUV uv;
			in >> uv.u() >> uv.v();
			if (in.fail())
			{
				std::cerr << "Mesh::read - Error reading uv coords.";
				clear();
				return false;
			}
			else if(uv.u()>1||uv.v()>1)
			{
				std::cerr << "Mesh::read - Error uv coords out of range";
				clear();
				return false;
			}
			tmp.push_back(uv);
		}
		m_textureArrays.push_back(tmp);
	}

	in >> str;
	if (str.compare("INDEXARRAY") != 0)
	{
		std::cerr << "Mesh::read - Expected INDEXARRAY directive found " << str;
		clear();
		return false;
	}

	m_indexArray.reserve(indices);
	for(;indices>0;indices--)
	{
		IndexedTri tri;

		in >> tri.a() >> tri.b() >> tri.c();

		if (in.fail())
		{
			std::cerr << "Mesh::read - Error reading indices";
			clear();
			return false;
		}
		m_indexArray.push_back(tri);
	}

	in >> str >> i;
	if (in.fail() || str.compare("FRAMES") != 0)
	{
		std::cerr << "Mesh::read - Expected FRAMES directive found " << str;
		clear();
		return false;
	}

	if (i > 0)
	{
		// TODO: animation frames
		for(; i > 0; --i)
		{
			in >> f >> f >> f >> f >> f >> f >> f >> f;
		}
		if (in.fail())
		{
			std::cerr << "Mesh::read - Error reading frames";
			clear();
			return false;
		}
	}

	in >> str >> i;
	if (in.fail() || str.compare("CONNECTORS") != 0)
	{
		std::cerr << "Mesh::read - Expected CONNECTORS directive found " << str;
		clear();
		return false;
	}

	if ( i > 0)
	{
		// TODO: Connectors
		for(; i > 0; --i)
		{
			in >> str >> f >> f >> f >> f >> f >> f ;
		}
		if(in.fail())
		{
			std::cerr << "Mesh::read - Error reading connectors";
			clear();
			return false;
		}
	}

	recalculateBoundData();

	return true;
}

void Mesh::write(std::ostream &out) const
{
	out << "MESH ";
	if (m_name.empty())
	{
		out << "_noname_\n";
	}
	else
	{
		out << m_name << '\n';
	}

	// noboolalpha should be default...
	out << "TEAMCOLOURS " << std::noboolalpha << teamColours() << '\n';

	out << "VERTICES " << vertices() << '\n';
	out << "FACES " << faces() << '\n';
	out << "VERTEXARRAY\n" ;

	std::vector<WZMVertex>::const_iterator vertIt;
	for (vertIt=m_vertexArray.begin(); vertIt < m_vertexArray.end(); vertIt++ )
	{
		out << '\t';
		out		<< vertIt->x() << ' '
				<< vertIt->y() << ' '
				<< vertIt->z() << '\n';
	}

	out << "TEXTUREARRAYS " << textureArrays() << '\n';

	std::vector< std::vector<WZMUV> >::const_iterator texArrIt;
	for (texArrIt=m_textureArrays.begin(); texArrIt < m_textureArrays.end(); texArrIt++ )
	{
		out << "TEXTUREARRAY " << std::distance(m_textureArrays.begin(),texArrIt) << '\n';
		std::vector<WZMUV>::const_iterator texIt;
		for (texIt=texArrIt->begin(); texIt < texArrIt->end(); texIt++ )
		{
			out << '\t';
			out	<< texIt->u() << ' '
					<< texIt->v() << '\n';
		}
	}

	out << "INDEXARRAY\n";

	std::vector<IndexedTri>::const_iterator indIt;
	for (indIt=m_indexArray.begin(); indIt < m_indexArray.end(); indIt++ )
	{
		out << '\t';
		out		<< indIt->a() << ' '
				<< indIt->b() << ' '
				<< indIt->c() << '\n';
	}

	// TODO: Frames and connectors
	out <<"FRAMES 0\nCONNECTORS 0\n";
}

bool Mesh::importFromOBJ(const std::vector<OBJTri>&	faces,
						 const std::vector<OBJVertex>& verts,
						 const std::vector<OBJUV>&	uvArray)
{
	const WZMVertex::less_wEps vertLess; // For make_mypair
	const WZMUV::less_wEps uvLess;

	typedef std::set<mypair<WZMVertex, WZMUV,
		WZMVertex::less_wEps, WZMUV::less_wEps> > t_pairSet;

	t_pairSet pairSet;

	std::vector<OBJTri>::const_iterator itF;

	std::pair<t_pairSet::iterator, bool> inResult;

	std::vector<unsigned> mapping;
	std::vector<unsigned>::iterator itMap;

	unsigned i;

	IndexedTri tmpTri;
	WZMUV tmpUv;

	clear();

	m_textureArrays.push_back(TexArray());

	for (itF = faces.begin(); itF != faces.end(); ++itF)
	{
		for (i = 0; i < 3; ++i)
		{
			/* in the uv's, -1 is "not specified," but the OBJ indices
			 * are 0 based, hence < 1
			 */
			if (itF->uvs.operator [](i) < 1)
			{
				tmpUv.u() = 0;
				tmpUv.v() = 0;
			}
			else
			{
				tmpUv = uvArray[itF->uvs.operator [](i)-1];
			}
			inResult = pairSet.insert(make_mypair(verts[itF->tri[i]-1],
												  tmpUv,
												  vertLess,
												  uvLess));

			if (!inResult.second)
			{
				tmpTri[i] = mapping[std::distance(pairSet.begin(), inResult.first)];
			}
			else
			{
				itMap = mapping.begin();
				std::advance(itMap, std::distance(pairSet.begin(), inResult.first));
				mapping.insert(itMap, m_vertexArray.size());
				tmpTri[i] = m_vertexArray.size();
				m_vertexArray.push_back(verts[itF->tri[i]-1]);
				m_textureArrays[0].push_back(tmpUv);
			}
		}
		m_indexArray.push_back(tmpTri);
	}

	recalculateBoundData();

	return true;
}

std::stringstream* Mesh::exportToOBJ(const Mesh_exportToOBJ_InOutParams& params) const
{
	const bool invertV = true;
	std::stringstream* out = new std::stringstream;

	std::pair<std::set<OBJVertex, OBJVertex::less_wEps>::iterator, bool> vertInResult;
	std::pair<std::set<OBJUV, OBJUV::less_wEps>::iterator, bool> uvInResult;

	std::vector<IndexedTri>::const_iterator itF;
	std::vector<unsigned>::iterator itMap;
	unsigned i;

	OBJUV uv;

	*out << "o " << m_name << "\n\n";

	for (itF = m_indexArray.begin(); itF != m_indexArray.end(); ++itF)
	{
		*out << "f";

		for (i = 0; i < 3; ++i)
		{
			*out << ' ';

			vertInResult = params.vertSet->insert(m_vertexArray[itF->operator [](i)]);

			if (!vertInResult.second)
			{
				*out << (*params.vertMapping)[std::distance(params.vertSet->begin(), vertInResult.first)] + 1;
			}
			else
			{
				itMap = params.vertMapping->begin();
				std::advance(itMap, std::distance(params.vertSet->begin(), vertInResult.first));
				params.vertMapping->insert(itMap, params.vertices->size());
				params.vertices->push_back(m_vertexArray[itF->operator [](i)]);
				*out << params.vertices->size();
			}

			*out << '/';

			uv = m_textureArrays[0][itF->operator [](i)];
			if (invertV)
			{
				uv.v() = 1 - uv.v();
			}
			uvInResult = params.uvSet->insert(uv);

			if (!uvInResult.second)
			{
				*out << (*params.uvMapping)[std::distance(params.uvSet->begin(), uvInResult.first)] + 1;
			}
			else
			{
				itMap = params.uvMapping->begin();
				std::advance(itMap, std::distance(params.uvSet->begin(), uvInResult.first));
				params.uvMapping->insert(itMap, params.uvs->size());
				params.uvs->push_back(uv);
				*out << params.uvs->size();
			}
		}
		*out << '\n';
	}

	return out;
}

std::string Mesh::getName() const
{
	return m_name;
}

void Mesh::setName(const std::string& name)
{
	if (isValidWzName(name))
	{
		m_name = name;
	}
}

bool Mesh::teamColours() const
{
	return m_teamColours;
}

void Mesh::setTeamColours(bool tc)
{
	m_teamColours = tc;
}

const WZMConnector& Mesh::getConnector(int index) const
{
	std::list<WZMConnector>::const_iterator pos;
	std::advance(pos, index);
	return *pos;
}

void Mesh::addConnector (const WZMConnector& conn)
{
	m_connectors.push_back(conn);
}

void Mesh::rmConnector (int index)
{
	int i;
	std::list<WZMConnector>::iterator pos;
	for(i=0,pos=m_connectors.begin();i<index;i++,pos++);
	if(pos==m_connectors.end())
	{
		return;
	}
	m_connectors.erase(pos);
}

int Mesh::connectors() const
{
	return m_connectors.size();
}

int Mesh::textureArrays() const
{
	return m_textureArrays.size();
}

const TexArray& Mesh::getTexArray (int index) const
{
	return m_textureArrays.at(index);
}

#if 0
void Mesh::addTexArray (const TexArray& tex, int index)
{
	if(tex.size()!=indices())
	{
		return;
	}
	m_textureArrays.insert(m_textureArrays.begin() + index,tex);
}
#endif
#if 0
void Mesh::rmTexArray(int index)
{
	std::vector<TexArray>::iterator pos;
	pos=m_textureArrays.begin()+index;
	if(pos==m_textureArrays.end())
	{
		return;
	}
	m_textureArrays.erase(pos);
}
#endif

unsigned Mesh::vertices() const
{
	return m_vertexArray.size();
}

unsigned Mesh::faces() const
{
	return triangles();
}

unsigned Mesh::triangles() const
{
	return m_indexArray.size();
}

unsigned Mesh::frames() const
{
	return m_frameArray.size();
}

unsigned Mesh::indices() const
{
	return m_indexArray.size();
}

bool Mesh::isValid() const
{
	// TODO: check m_frameArray, m_connectors
	if (!isValidWzName(m_name))
	{
		return false;
	}

	// Check that the values of the indices are in range
	std::vector<IndexedTri>::const_iterator it;
	for (it = m_indexArray.begin(); it != m_indexArray.end(); ++it)
	{
		if ((*it).a() >= m_vertexArray.size())
		{
			return false;
		}
		if ((*it).b() >= m_vertexArray.size())
		{
			return false;
		}
		if ((*it).c() >= m_vertexArray.size())
		{
			return false;
		}
	}
	return true;
}

void Mesh::defaultConstructor()
{
	m_name.clear();
	m_teamColours = false;
}

void Mesh::clear()
{
	m_name.clear();
	m_frameArray.clear();
	m_vertexArray.clear();
	m_textureArrays.clear();
	m_indexArray.clear();
	m_connectors.clear();
	m_teamColours = false;
}

void Mesh::scale(GLfloat x, GLfloat y, GLfloat z)
{
	std::vector<WZMVertex>::iterator vertIt;
	for (vertIt = m_vertexArray.begin(); vertIt < m_vertexArray.end(); ++vertIt )
	{
		vertIt->scale(x, y, z);
	}

	std::list<WZMConnector>::iterator itC;
	for (itC = m_connectors.begin(); itC != m_connectors.end(); ++itC)
	{
		itC->m_pos.scale(x, y, z);
	}

	m_mesh_weightcenter.scale(x, y, z);
	m_mesh_aabb_min.scale(x, y, z);
	m_mesh_aabb_max.scale(x, y, z);

}

void Mesh::mirrorUsingLocalCenter(int axis)
{
	mirrorFromPoint(getCenterPoint(), axis);
}

void Mesh::mirrorFromPoint(const WZMVertex& point, int axis)
{
	std::vector<WZMVertex>::iterator vertIt;
	for (vertIt = m_vertexArray.begin(); vertIt < m_vertexArray.end(); ++vertIt )
	{
		switch (axis)
		{
		case 0:
			vertIt->x() = -vertIt->x() + 2 * point.x();
			break;
		case 1:
			vertIt->y() = -vertIt->y() + 2 * point.y();
			break;
		default:
			vertIt->z() = -vertIt->z() + 2 * point.z();
		}
	}

	std::list<WZMConnector>::iterator itC;
	for (itC = m_connectors.begin(); itC != m_connectors.end(); ++itC)
	{
		switch (axis)
		{
		case 0:
			itC->m_pos.operator[](0) = -itC->m_pos.operator[](0) + 2 * point.x();
			break;
		case 1:
			itC->m_pos.operator[](1) = -itC->m_pos.operator[](1) + 2 * point.y();
			break;
		default:
			itC->m_pos.operator[](2) = -itC->m_pos.operator[](2) + 2 * point.z();
		}
	}

	// for convenience
	reverseWinding();

	recalculateBoundData();
}

void Mesh::reverseWinding()
{
	std::vector<IndexedTri>::iterator it;
	for (it = m_indexArray.begin(); it != m_indexArray.end(); ++it)
	{
		std::swap((*it).b(), (*it).c());
	}
}

void Mesh::recalculateBoundData()
{
	WZMVertex weight, min, max;

	if (m_vertexArray.size())
	{
		min = max = m_vertexArray.at(0);

		std::vector<WZMVertex>::const_iterator vertIt;
		for (vertIt = m_vertexArray.begin(); vertIt < m_vertexArray.end(); ++vertIt )
		{
			weight.x() += vertIt->x();
			weight.y() += vertIt->y();
			weight.z() += vertIt->z();

			if (min.x() > vertIt->x()) min.x() = vertIt->x();
			if (min.y() > vertIt->y()) min.y() = vertIt->y();
			if (min.z() > vertIt->z()) min.z() = vertIt->z();

			if (max.x() < vertIt->x()) max.x() = vertIt->x();
			if (max.y() < vertIt->y()) max.y() = vertIt->y();
			if (max.z() < vertIt->z()) max.z() = vertIt->z();
		}

		weight.x() /= m_vertexArray.size();
		weight.y() /= m_vertexArray.size();
		weight.z() /= m_vertexArray.size();
	}

	m_mesh_weightcenter = weight;
	m_mesh_aabb_min = min;
	m_mesh_aabb_max = max;
}

WZMVertex Mesh::getCenterPoint() const
{
	WZMVertex center;

	center.x() = (m_mesh_aabb_max.x() + m_mesh_aabb_min.x()) / 2;
	center.y() = (m_mesh_aabb_max.y() + m_mesh_aabb_min.y()) / 2;
	center.z() = (m_mesh_aabb_max.z() + m_mesh_aabb_min.z()) / 2;

	return center;
}
