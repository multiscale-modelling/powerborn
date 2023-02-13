/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE for License
 *    and Copyright information.
 */

#ifndef POWERSASA2_H_
#define POWERSASA2_H_

#include "Typedefs.h"
#include "PowerCell.h"

#include <cstdlib>
#include <cmath>
#include <iostream>

#define MAX_NB 20
#define MAX_VX 12
#define MAX_PNT 4
#define MAX_COUNT 100

namespace powerborn
{

class PowerSasa2Exception: public std::exception
{
};

typedef Coord3 Coord;
typedef AlignedVector<Coord>::Type CoordVec;
typedef AlignedVector<float>::Type FloatVec;
typedef AlignedVector<int>::Type IntVec;

class PowerSasa2
{

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    inline float DRAD2()
    {
        return static_cast<float> (1000.0)
                * std::numeric_limits<float>::epsilon();
    }
    inline float DANG()
    {
        return static_cast<float> (1000.0)
                * std::numeric_limits<float>::epsilon();
    }
    inline float pi()
    {
        return static_cast<float> (3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342);
    } //everything const, will be optimized away.
    inline float calc_sasa_single(PowerCell& atom);

    PowerSasa2(): do_debug_(0), debug_aid_(0)
    {
        Init();
    }
    PowerSasa2(unsigned int aid): do_debug_(1), debug_aid_(aid)
    {
        std::cout << "PowerSasa2 debugging atom " << aid << std::endl;
        Init();
    }
private:
    inline void Init()
    {
        Resize_NB(MAX_NB);
        Resize_VX(MAX_VX);
        Resize_PNT(MAX_PNT);
        Resize_NA();
    }
    inline void Resize_NB(unsigned int nnb)
    {
        np.resize(nnb);
        e.resize(nnb);
        sintheta.resize(nnb);
        costheta.resize(nnb);
        unsigned int npnt = MAX_PNT;
        if (p.size() > 0)
            npnt = p[0].size();
        unsigned int nnb_old = p.size();
        next.resize(nnb);
        p.resize(nnb);
        ang.resize(nnb);
        for (unsigned int i = nnb_old; i < nnb; ++i)
        {
            next[i].resize(npnt);
            p[i].resize(npnt);
            ang[i].resize(npnt);
        }
        //std::cout << "db> resize nnb=" << nnb << std::endl;
    }
    inline void Resize_VX(unsigned int nvx)
    {
        off.resize(nvx);
        vx.resize(nvx);
        br_c.resize(nvx);
        br_p.resize(nvx);
        //std::cout << "db> resize nvx=" << nvx << std::endl;
    }
    inline void Resize_PNT(unsigned int npnt)
    {
        for (unsigned int i = 0; i < p.size(); ++i)
        {
            next[i].resize(npnt);
            p[i].resize(npnt);
            ang[i].resize(npnt);
        }
        rang.resize(npnt);
        pos.resize(npnt);
        //std::cout << "db> resize npnt=" << npnt << std::endl;
    }
    inline void Resize_NA()
    {
        tol_pow = DRAD2();
    }

    inline float Ang_About(Coord const& a, Coord const& b, Coord const& c);
    inline void Get_Ang(const int & np, const IntVec & p, const Coord & e,
            const float & sintheta, const float & costheta,
            FloatVec & ang);
    inline void Get_Next(int n, FloatVec & ang, IntVec & next,
            const IntVec & p, const Coord & e);
    inline int get_Neighbour_Properties(const PowerCell& pc);
    inline int register_Surface_Vertices(const PowerCell& pc);
    inline float sasa_from_contours(const int nvx, unsigned int aid);
    inline void sasa_from_single_circles(float& sasa_ia, const PowerCell& pc);

    float tol_pow;

    // ----- for calc_sasa_single --------------------

    IntVec np; // number of points (registered vertices) of i-th atom
    CoordVec e; // direction to neighbour
    FloatVec sintheta;
    FloatVec costheta;

    IntVec off; // off[i] != 0 if vertix i is already taken into account
    CoordVec vx; // all surface vertices
    AlignedVector<UInt2>::Type br_c; // bridge between circles (intersections with neighbors)
    AlignedVector<UInt2>::Type br_p; // bridge between point numbers

    AlignedVector<FloatVec>::Type ang; // angles
    AlignedVector<IntVec>::Type next; // next[][i] is the number of angle that follows the i-th one
    AlignedVector<IntVec>::Type p; // number of surface vertices

    //----- for Get_Next ----------------------------

    IntVec rang;
    IntVec pos;

    const unsigned int do_debug_, debug_aid_;
};
//Implementation:

//-----------------------------------------------------------------------------
// Retuns twist angle (in interval [0, 2*pi()) ) between Vectors a and b about c.
// a and b should be of unit length and perpendicular to c.

float PowerSasa2::Ang_About(
        Coord const& a, Coord const& b, Coord const& c)
{
    float const ONE = 0.999f;
    float const THRESHOLD = 0.001f;
    float ang, co, vp;
    Coord v;

    co = a.dot(b);
    if (co <= -ONE)
    {
        v = a.cross(b);
        ang = pi() - std::asin(v.norm());
    }
    else if (co >= ONE)
    {
        v = a.cross(b);
        ang = std::asin(v.norm());
    }
    else
        ang = std::acos(co);

    if (std::fabs(c[0]) > THRESHOLD)
    {
        vp = a[1] * b[2] - a[2] * b[1];
        if ((vp < 0.0f) != (c[0] < 0.0f))
            ang = -ang;
    }
    else if (std::fabs(c[1]) > THRESHOLD)
    {
        vp = a[2] * b[0] - a[0] * b[2];
        if ((vp < 0.0f) != (c[1] < 0.0f))
            ang = -ang;
    }
    else if (std::fabs(c[2]) > THRESHOLD)
    {
        vp = a[0] * b[1] - a[1] * b[0];
        if ((vp < 0.0f) != (c[2] < 0.0f))
            ang = -ang;
    }
    else
    {
        std::cerr << "PowerSasa2: Axis too short" << std::endl;
        throw PowerSasa2Exception();
    }
    if (ang < 0.0f)
        ang += 2.0f * pi();
    return ang;
}
//-----------------------------------------------------------------------------
//Finds phi-angles of s-vertices

void PowerSasa2::Get_Ang(
        const int & np, // total number of s-vertices for given circle
        const IntVec & p, // s-vertice number for given circle
        const Coord & e, // direction to neighbour
        const float & sintheta, const float & costheta, // theta angle
        FloatVec & ang) // output

{
    //    if (np == 0) return;//is tested outside
    ang[0] = 0.0;
    float inv_sintheta = 1.0f / sintheta;
    Coord costheta_e = costheta * e;
    Coord pu0 = (vx[p[0]] - costheta_e) * inv_sintheta;
    for (int j = 1; j < np; j++)
    {
        Coord pu = (vx[p[j]] - costheta_e) * inv_sintheta;
        ang[j] = Ang_About(pu0, pu, e);
    }
}

//-----------------------------------------------------------------------------
// given array ang of length n founds number next[i] of ang element that is next in magnitude to element i
// can be optimized (???)

void PowerSasa2::Get_Next(
        int n, FloatVec & ang, IntVec & next,
        const IntVec & p, const Coord & e)
{
    int j, k, m;
    //static int rang[MAX_PNT], pos[MAX_PNT];
    if (n == 0)
        return;

    for (j = 1; j < n; ++j)
    {

        if (ang[j] <= 2.0f * pi() - DANG())
            continue;
        float dp = vx[p[j]].cross(vx[p[0]]).dot(e);
        if (dp < 0.0f)
            ang[j] = 0.0f;
        else if (dp == 0.0f)
        {
            std::cout << "PowerSasa2: Precision insufficient to resolve angles 1"
                    << std::endl;
            Eigen::IOFormat fmt(4, 0, "", "", " ", "", "(", ")");
            std::cout <<  p[j] << " | " << p[0] << " | " <<  j << " | " << 0 << std::endl;
            std::cout << vx[p[j]].format(fmt) << " X " << vx[p[0]].format(fmt) << " XX " << e.format(fmt) << " X " << std::endl;
            std::cout << (vx[p[j]] - vx[p[0]]).norm() << std::endl;
            throw PowerSasa2Exception();
        }
    }

    rang[0] = 0;
    rang[1] = 1;

    for (j = 2; j < n; j++)
    {
        m = j;
        for (k = 1; k < j; k++)
        {
            if (ang[k] > ang[j] + DANG())
            {
                (rang[k])++;
                m--;
            }
            else if (ang[k] > ang[j] - DANG())
            {

                float dp = vx[p[j]].cross(vx[p[k]]).dot(e);
                if (dp > 0.0f)
                {
                    (rang[k])++;
                    m--;
                }
                else if (dp == 0.0f)
                {
                    //TODO why is this???
                    std::cout
                            << "PowerSasa2: Precision insufficient to resolve angles 2"
                            << std::endl;
                    Eigen::IOFormat fmt(4, 0, "", "", " ", "", "(", ")");
                    std::cout <<  p[j] << " | " << p[k] << " | " <<  j << " | " << k << std::endl;
                    std::cout <<  vx[p[j]].format(fmt) << " X " << vx[p[k]].format(fmt) << " X " << e.format(fmt) << std::endl;
                    throw PowerSasa2Exception();
                }
            }
        }
        rang[j] = m;
    }

    for (j = 0; j < n; j++)
        pos[rang[j]] = j;

    for (j = 0; j < n; j++)
    {
        if (rang[j] == n - 1)
            next[j] = pos[0];
        else
            next[j] = pos[rang[j] + 1];
    }
}

//===================================================================================

int PowerSasa2::get_Neighbour_Properties(const PowerCell& pc)
{
    unsigned int n_apart = 0;
    float radius  = pc.normedRadius();
    float radius2 = radius * radius;
    float inv_rad = 1.0f / radius;
    unsigned int nnb = pc.neighbourSize();

    bool invalid_dist = 0;
    bool covered = 0;
    bool cos_theta_error = 0;
    for (unsigned int i = 0; i < nnb; ++i)       // over neighbors
    {
        Coord rel_pos = pc.position(i);       // vector to neighbor
        float dist = rel_pos.norm();                  // distance to neighbor
        invalid_dist |= dist == 0.0f;
        dist = dist == 0.0f ? 1e-6f : dist;
        float nb_rad = pc.radius(i);                     // neighbor RADius
        covered |= dist <= nb_rad - radius;
        float inv_dist = 1.0f / dist;
        costheta[i] = 0.5f * ( dist + (radius2 - nb_rad * nb_rad) * inv_dist ) * inv_rad;
        cos_theta_error |= costheta[i] <= -1.0f;
        np[i] = 0;
        sintheta[i] = std::sqrt(1.0f - costheta[i] * costheta[i]);
        e[i] = rel_pos * inv_dist;
    }
    if(invalid_dist)
    {
            std::cerr << "PowerSasa: Invalid distance to neighbour" << std::endl;
            throw PowerSasa2Exception();
    }
    bool rv = covered | cos_theta_error | (n_apart==nnb);
    return rv;
}

//===================================================================================

int PowerSasa2::register_Surface_Vertices(const PowerCell& pc)
{
    int nvx = pc.zeroPointSize();
    if (nvx >= static_cast<int> (vx.size()))
    {
        Resize_VX(nvx);
    }
    {
        for (unsigned int k = 0; k < pc.zeroPointSize(); ++k)
        {
            const ZeroPoint& zp = pc.zeroPoint(k);
            // the two common generators
            int ptn0 = zp.generators()[0];
            int ptn1 = zp.generators()[1];

            if ((np[ptn0] >= static_cast<int> (rang.size())) || (np[ptn1]
                                                                    >= static_cast<int> (rang.size())))
            {
                Resize_PNT((np[ptn0] > np[ptn1]) ? np[ptn0] + 1 : np[ptn1] + 1);
                //std::cerr << "PowerSasa2: Number of surface vertices for single neighbor exeeds MAX_PNT=" << MAX_PNT << std::endl;
                //throw PowerSasa2Exception();
            }

            vx[k] = zp.position();
            p[ptn0][np[ptn0]] = p[ptn1][np[ptn1]] = k; // s-vertex number

            br_c[k][0] = ptn0;
            br_c[k][1] = ptn1;
            br_p[k][0] = np[ptn0];
            br_p[k][1] = np[ptn1];
            ++np[ptn0];
            ++np[ptn1];
        }
    }
    return nvx;
}

//===================================================================================

float PowerSasa2::sasa_from_contours(const int nvx, unsigned int aid)
{
    float sasa_ia = 0.0f;
    unsigned int ic1, ic2, ic_0, ic_1, ip2, ip_next, ivx, count;
    float phi, co, dirdet;
    Coord *p_ini, *pt,  vv;

    for (int iv = 0; iv < nvx; iv++)
        off[iv] = 0;

    for (int iv = 0; iv < nvx; iv++)
    {
        if (off[iv])
            continue;
        p_ini = &vx[iv];
        ic_0 = br_c[iv][0];
        ic_1 = br_c[iv][1];

        // determine which directoin
        dirdet = (e[ic_1].cross(e[ic_0])).dot(*p_ini);
        if (dirdet == 0.0f)
        {
            std::cerr << "PowerSasa2: dirdet == 0.0" << std::endl;
            throw PowerSasa2Exception();
        }
        //assert(dirdet != 0.0);
        bool dirdet_gt0 = dirdet > 0.0f;
        ic1 = dirdet_gt0? ic_0 : ic_1;
        ic2 = dirdet_gt0? ic_1 : ic_0;
        ip2 = dirdet_gt0? br_p[iv][1] : br_p[iv][0];

        /*if(do_debug_ && aid == debug_aid_)
        {
            std::cout << " NEW dirdet " << dirdet << " " << " " << std::endl;
        }*/
        pt = p_ini;
        sasa_ia += 2.0f * pi();
        count = 0;

        do // main loop
        {
            ++count;
            if (count > MAX_COUNT)
            {
                std::cerr << "PowerSasa2: Wrong contour" << std::endl;
                throw PowerSasa2Exception();
            }
            //assert(count <= MAX_PNT);   // else wrong contour
            ip_next = next[ic2][ip2];
            phi = ang[ic2][ip_next] - ang[ic2][ip2];
            phi += phi < 0.0f ? (2.0f * pi()) : 0.0f;
            //if (phi < 0.0f) phi += 2.0f * pi();

            co = (e[ic1].dot(e[ic2]) - costheta[ic1] * costheta[ic2])
                                        / (sintheta[ic1] * sintheta[ic2]);

            co = co < -1.0f ? -1.0f : co;
            co = co >  1.0f ?  1.0f : co;
            float t_sasa = phi * costheta[ic2] - acos(co);
            if(t_sasa > 0.0f && (fabs(phi - 2.0f * pi()) < 1e-4f)) // numerical stability
            {
                //std::cout << "PowerSasa2 Warning: very small phi angle !!!" << std::endl;
                phi = -(phi - 2.0f * pi());
                t_sasa = phi * costheta[ic2] - acos(co);
            }
            sasa_ia += t_sasa;

            if(do_debug_ && debug_aid_ == aid)
                std::cout << "NEW " << sasa_ia << " t_sasa " << t_sasa
                    << " phi " << phi << " co " << co << " costheta ic2 " << costheta[ic2]
                    << " angic2ip2 " << ang[ic2][ip2] << " angic2next " << ang[ic2][ip_next]
                    << std::endl;

            off[p[ic2][ip2]] = 1;
            ic1 = ic2;

            ivx = p[ic1][ip_next];
            pt = &vx[ivx];
            bool tmp = br_c[ivx][0] == ic1;
            ic2 = tmp ? br_c[ivx][1] : br_c[ivx][0];
            ip2 = tmp ? br_p[ivx][1] : br_p[ivx][0];

            /*if (br_c[ivx][0] == ic1)
            {
                ic2 = br_c[ivx][1];
                ip2 = br_p[ivx][1];
            }
            else
            {
                ic2 = br_c[ivx][0];
                ip2 = br_p[ivx][0];
            }*/
        } while (pt != p_ini);
        if(do_debug_ && debug_aid_ == aid)
        {
            std::cout << "NEW Contour with " << count << " surface vertices has sasa "
                    << sasa_ia / pi() / 4.0f << std::endl;
        }
        if (sasa_ia > 4.0f * pi())
            sasa_ia -= 4.0f * pi();
    }
    return sasa_ia;
}

//===================================================================================

void PowerSasa2::sasa_from_single_circles(float& sasa_ia, const PowerCell& pc)
{
    int nnb = pc.neighbourSize();
    for (int i = 6; i < nnb; ++i) // skip the first six dummy neighbours
    {
        if ((np[i] != 0) | pc.singleCircleExcluded(i))
            continue;
        float dist = pc.position(i).norm();

        if (dist >= pc.normedRadius() + pc.radius(i) || dist <= pc.normedRadius() - pc.radius(i))
            continue;

        int ok = 1;
        Coord cc = costheta[i] * e[i] * pc.normedRadius();
        //float pw_i = -sintheta[i] * sintheta[i]; // this is numerically unstable due to the std::sqrt used to compute sintheta
        float pw_i = cc.squaredNorm() - pc.normedRadius() * pc.normedRadius(); // so just recompute the correct value
        for (int j = 6; j < nnb; ++j) // skip the first six dummy neighbours
        {
            float pw_j = (pc.position(j) - cc).squaredNorm() - pc.radius(j) * pc.radius(j);
            if (pw_j <= pw_i && i != j)
            {
                ok = 0;
                break;
            }
        }
        if (ok)
        {
            sasa_ia += 2.0f * pi() * (1.0f + costheta[i]);
            if (sasa_ia > 4.0f * pi())
                sasa_ia -= 4.0f * pi();
        }
    }
}

//===================================================================================

float PowerSasa2::calc_sasa_single(PowerCell& pc)
{
    const int nnb = pc.neighbourSize(); // number of neighbors
    //std::cout << "Sasa id " << pc.id() << std::endl;
    if(pc.isCovered())
    {
        //if(do_debug_ && pc.id() == debug_aid_) std::cout << "PowerSasa is covered" << std::endl;
        return 0.0f;
    }
    else if(nnb <= 6) // only dummy neighbours
    {
        //if(do_debug_ && pc.id() == debug_aid_) std::cout << "PowerSasa no neighbours" << std::endl;
        return 4.0f * pi() * pc.radius() * pc.radius() * (pc.normedRadius() * pc.normedRadius());
    }

    if (nnb >= static_cast<int> (np.size()))
    {
        Resize_NB(nnb + 1);
    }

    //----- get angles and other properties of neighbors -----

    int do_return = this->get_Neighbour_Properties(pc);
    if(do_return)
    {
        //if(do_debug_ && pc.id() == debug_aid_) std::cout << "PowerSasa no overlapping neighbours" << std::endl;
        if(nnb > 6)
        {
            // in Diagram shpere has neighbours that could cut it, but here no neighbours actually cut
            // so it must be completely inside another sphere!
            return 0.0f;
        }
        return 4.0 * pi() * pc.radius() * pc.radius() * (pc.normedRadius() * pc.normedRadius());
    }

    //------ register surface vertices ------------------------
    int nvx = this->register_Surface_Vertices(pc);

    // ----- get angles of s-vertices within each circle ---------
    for (int i = 0; i < nnb; i++)
    {
        if (np[i] <= 0)
            continue;
        if (np[i] % 2 != 0)
        {
            std::cerr << "PowerSasa2: odd number of crossings for Cell " << pc.id()
                << " nb_id " << i << " zeropoints " << np[i] << std::endl;
            //pc.writeZeroPoints("odd_crossing_zeropoints.pqr");
            throw PowerSasa2Exception();
        }
        //if(do_debug_ && pc.id() == debug_aid_) std::cout << "circle " << i << std::endl;
        Get_Ang(np[i], p[i], e[i], sintheta[i], costheta[i], ang[i]);
        Get_Next(np[i], ang[i], next[i], p[i], e[i]);
    }

    // ---------- get sasa from contours ---------------------------
    float sasa_ia = this->sasa_from_contours(nvx, pc.id());
    //if(do_debug_ && pc.id() == debug_aid_)
    //    std::cout << "NEW sasa from contours " << sasa_ia / pi() / 4.0f<< std::endl;
    // ---------- sasa from single circles ----------------------------
    this->sasa_from_single_circles(sasa_ia, pc);
    //if(do_debug_ && pc.id() == debug_aid_)
    //   std::cout << "NEW sasa with single circles " << sasa_ia / pi() / 4.0f << std::endl;
    sasa_ia = sasa_ia < 0.0f ? 0.0f : sasa_ia; // catch numerical problems
    return sasa_ia * pc.radius() * pc.radius() * (pc.normedRadius() * pc.normedRadius());
}

//=============================================================

}

#endif /* POWERSASA2_H_ */
