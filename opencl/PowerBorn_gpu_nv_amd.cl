/*
 *    This file is part of the PowerBorn Library
 *    Please refer to the file LICENSE (in the root of the 
 *    PowerBorn include directory) for License
 *    and Copyright information.
 */

#ifndef FIT_PARAM_A
#define FIT_PARAM_A 1.0f
#endif

#ifndef FIT_PARAM_B
#define FIT_PARAM_B 0.0f
#endif

#ifndef GPU_BLOCK_SIZE
#define GPU_BLOCK_SIZE 256
#endif

#ifndef GPU_BLOCK_SIZE_KERNEL2
#define GPU_BLOCK_SIZE_KERNEL2 128
#endif

#define M_PI 3.141592653589793f


#define MAX_GPU_MAINCELLCOUNT 512
#define LAST_LEVEL_GRID_SIZE 32
#define WARP_SIZEL 64
#define WARP_SIZE 32
#define WARP_SIZE2 16
#define LEVEL3_SIZE 64
#define LEVEL3_SHIFT 6
#define CUR_ATOMS_SIZE 32
#define SQRT3 1.732052f
#define WATERMASK_SIZE 512

#define WATERMASK_SMALL 256
#define WATERMASK_LARGE 512
#define WATERMASK_THRESHOLD 3.00001f

typedef struct
{
    float4 root;
    unsigned int atomCnt;
    unsigned int atomPitch;
    unsigned int waterCnt;
    unsigned int atomsOffset;
    unsigned int waterOffset;
    unsigned int stackOffset;
    unsigned int stackSize;
    //unsigned int waterSeedsCnt;
    unsigned int octTreeMaxSize;
    float integration_factor;
    unsigned int maxLevel;
    unsigned int maincelllevel;
    unsigned int maincellcount; // count of non-zero main cells. this is the ammount of octtree nodes to start integrating with
} pGPUConstants;

typedef struct
{
    unsigned int dihCnt;
    float coulomb_param;
    float gb_param;
    float npse_param;
    float excluded_param_co;
    float excluded_param_lj;
    float r1, r2;
} pGPUEnergyConstants;

typedef struct
{
    float4 p;
    float r1,r2,r3,r4; // reserved to make the size == 96b
    int celltype[72];
} pPGUOctNode64;


float invR(float x, float y, float z)
{
    x*=x;
    y*=y;
    z*=z;
    return 1.0f / sqrt(x+y+z);
}

float squareIntegral(float x, float y, float z, float angley, float anglez)
{
    // used only for bounding box
    float xinv = 1.0f/x;
    float yt = y*xinv;
    float zt = z*xinv;
    float ytsq = yt*yt;
    float ztsq = zt*zt;
    float f1=1.0f+ytsq;
    float f2=1.0f+ztsq;
    float g1a=f1+f2;
    float g1=yt*zt;
    g1*=g1a;
    float g2a=f2+ytsq;
    float g2=f1*f2*g2a;
    float tmp1=sqrt(f1);
    float test3=anglez;
    test3 = zt < 0 ? -test3 : test3;
    float gb_1a=2.0f*yt*test3;
    float gb_1b=0.5f+f1;
    float gb_1=gb_1a*gb_1b;
    float gb_2=f1*tmp1;
    float tmp2=sqrt(f2);
    float test4=angley;
    test4 = yt < 0 ? -test4 : test4;
    float gc_1a=2.0f*zt*test4;
    float gc_1b=0.5f+f2;
    float gc_1=gc_1a*gc_1b;
    float gc_2=f2*tmp2;
    float gcb2 =gb_2*gc_2;
    float g_nom2=gb_1*g2*gc_2;
    float g_nom3=gc_1*g2*gb_2;
    g_nom2+=g_nom3;
    float g_nom1=g1*gcb2;
    float g_denom=g2*gcb2;
    g_nom2+=g_nom1;
    g_nom2 /= g_denom;
    float res2=xinv*xinv;
    const float res3=(-1.0f/24.0f)*xinv;
    return g_nom2*res2*res3;
}

void AssignTopLevelAtomAndWaterToCells(
        global int* stack,
        local int* render,
        const int bWrite,
        const float4 level1root,
        const int level1cellsperdimmask,
        const int level1cellsperdim,
        const int level1cellsperdim2,
        const int threadIdx,
        const pGPUConstants localConstants,
        global const float4* atoms,
        global const float4* water
)
{
    for(int atomId = threadIdx ; atomId < localConstants.atomCnt; atomId += GPU_BLOCK_SIZE)
    {
        float4 a = atoms[atomId + localConstants.atomsOffset];
        float4 d;
        d.xyz = a.xyz - level1root.xyz;
        d.w = a.w + level1root.w*SQRT3*0.5f;
        int px = d.x/level1root.w;
        int py = d.y/level1root.w;
        int pz = d.z/level1root.w;
        d.x -= px*level1root.w;
        d.y -= py*level1root.w;
        d.z -= pz*level1root.w;
        d.w *=d.w;
        // assume grid cell size is larger than radius
        int xxp = (px < level1cellsperdimmask && level1root.w < d.x + a.w);
        int xxm = -(px > 0 && a.w > d.x);
        int yyp = (py < level1cellsperdimmask && level1root.w < d.y + a.w);
        int yym = -(py > 0 && a.w > d.y);
        int zzp = (pz < level1cellsperdimmask && level1root.w < d.z + a.w);
        int zz = -(pz > 0 && a.w > d.z);
        int oz =  px+py*level1cellsperdim + (pz + zz)*level1cellsperdim2;
        d.xyz -= 0.5f*level1root.w;
        for(; zz <= zzp; zz++, oz+=level1cellsperdim2)
        {
            for(int yy = yym; yy <= yyp; yy++)
            {
                int o = xxm + yy*level1cellsperdim + oz;
                for(int xx = xxm; xx <= xxp; xx++, o++)
                {
                    if((d.x-xx*level1root.w)*(d.x-xx*level1root.w) +
                            (d.y-yy*level1root.w)*(d.y-yy*level1root.w) + 
                            (d.z-zz*level1root.w)*(d.z-zz*level1root.w) < d.w
                    )
                    {
                        unsigned int position = atomic_add(&render[o], 1) ;
                        if(bWrite)
                            stack[position + localConstants.stackOffset] = atomId;
                    }
                }
            }
        }
    }

    // assign water atoms to cells
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int atomId = threadIdx ; atomId < localConstants.waterCnt; atomId += GPU_BLOCK_SIZE)
    {
        float4 a = water[atomId + localConstants.waterOffset];
        float4 d;
        d.xyz = a.xyz - level1root.xyz;
        d.w = a.w + level1root.w*SQRT3*0.5f;
        int px = d.x/level1root.w;
        int py = d.y/level1root.w;
        int pz = d.z/level1root.w;
        d.x -= px*level1root.w;
        d.y -= py*level1root.w;
        d.z -= pz*level1root.w;
        d.w *=d.w;
        // assume grid size is larger than radius
        int xxp = (px < level1cellsperdimmask && level1root.w < d.x + a.w);
        int xxm = -(px > 0 && a.w > d.x);
        int yyp = (py < level1cellsperdimmask && level1root.w < d.y + a.w);
        int yym = -(py > 0 && a.w > d.y);
        int zzp = (pz < level1cellsperdimmask && level1root.w < d.z + a.w);
        int zz = -(pz > 0 && a.w > d.z);
        int oz =  px+py*level1cellsperdim + (pz + zz)*level1cellsperdim2;
        d.xyz -= 0.5f*level1root.w;
        for(; zz <= zzp; zz++, oz+=level1cellsperdim2)
        {
            for(int yy = yym; yy <= yyp; yy++)
            {
                int o = xxm + yy*level1cellsperdim + oz;
                for(int xx = xxm; xx <= xxp; xx++, o++)
                {
                    if(render[o] != 0)
                        if((d.x-xx*level1root.w)*(d.x-xx*level1root.w) +
                                (d.y-yy*level1root.w)*(d.y-yy*level1root.w) +
                                (d.z-zz*level1root.w)*(d.z-zz*level1root.w) < d.w
                        )

                        {
                            unsigned int position = atomic_add(&render[o + MAX_GPU_MAINCELLCOUNT], 1);
                            if(bWrite)
                                stack[position + localConstants.stackOffset] = atomId;
                        }
                }
            }
        }
    }
}

void AssignTopLevelAtomAndWaterToCellsLargeGrid(
        global int* stack,
        local int* render,
        const int bWrite,
        const float4 levelroot,
        const int levelcellsperdimmask,
        const int levelcellsperdim,
        const int levelcellsperdim2,
        const int threadIdx,
        const pGPUConstants localConstants,
        global const float4* atoms,
        global const float4* water
)
{
    for(int atomId = threadIdx ; atomId < localConstants.atomCnt; atomId += GPU_BLOCK_SIZE)
    {
        float4 a = atoms[atomId + localConstants.atomsOffset];
        float4 d;
        d.xyz = a.xyz - levelroot.xyz;
        d.w = a.w + levelroot.w*SQRT3*0.5f;
        int px = d.x/levelroot.w;
        int py = d.y/levelroot.w;
        int pz = d.z/levelroot.w;
        d.x -= px*levelroot.w;
        d.y -= py*levelroot.w;
        d.z -= pz*levelroot.w;
        d.w *=d.w;
        // assume grid size is larger than radius
        int xxp = (px < levelcellsperdimmask && levelroot.w < d.x + a.w);
        int xxm = -(px > 0 && a.w > d.x);
        int yyp = (py < levelcellsperdimmask && levelroot.w < d.y + a.w);
        int yym = -(py > 0 && a.w > d.y);
        int zzp = (pz < levelcellsperdimmask && levelroot.w < d.z + a.w);
        int zz = -(pz > 0 && a.w > d.z);
        d.xyz -= 0.5f*levelroot.w;
        for(; zz <= zzp; zz++)
        {
            for(int yy = yym; yy <= yyp; yy++)
            {
                for(int xx = xxm; xx <= xxp; xx++)
                {
                    if((d.x-xx*levelroot.w)*(d.x-xx*levelroot.w) +
                            (d.y-yy*levelroot.w)*(d.y-yy*levelroot.w) + 
                            (d.z-zz*levelroot.w)*(d.z-zz*levelroot.w) < d.w
                    )
                    {
                        int oLevel1 = ((px+xx) >> 2) + ((py+yy) >> 2)*(levelcellsperdim) + ((pz+zz) >> 2)*(levelcellsperdim2);
                        int o = ((px+xx) & 0x3) + (((py+yy) & 0x3)<< 2) + (((pz+zz) & 0x3) << 4);
                        unsigned int position = atomic_add(&render[oLevel1], 1) ;
                        if(bWrite)
                            stack[position + localConstants.stackOffset] = atomId + (o << 24);
                    }
                }
            }
        }
    }

    // assign water atoms to cells
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int atomId = threadIdx ; atomId < localConstants.waterCnt; atomId += GPU_BLOCK_SIZE)
    {
        float4 a = water[atomId + localConstants.waterOffset];
        float4 d;
        d.xyz = a.xyz - levelroot.xyz;
        d.w = a.w + levelroot.w*SQRT3*0.5f;
        int px = d.x/levelroot.w;
        int py = d.y/levelroot.w;
        int pz = d.z/levelroot.w;
        d.x -= px*levelroot.w;
        d.y -= py*levelroot.w;
        d.z -= pz*levelroot.w;
        d.w *=d.w;
        // assume grid size is larger than radius
        int xxp = (px < levelcellsperdimmask && levelroot.w < d.x + a.w);
        int xxm = -(px > 0 && a.w > d.x);
        int yyp = (py < levelcellsperdimmask && levelroot.w < d.y + a.w);
        int yym = -(py > 0 && a.w > d.y);
        int zzp = (pz < levelcellsperdimmask && levelroot.w < d.z + a.w);
        int zz = -(pz > 0 && a.w > d.z);
        d.xyz -= 0.5f*levelroot.w;
        for(; zz <= zzp; zz++)
        {
            for(int yy = yym; yy <= yyp; yy++)
            {
                for(int xx = xxm; xx <= xxp; xx++)
                {
                    int oLevel1 = ((px+xx) >> 2) + ((py+yy) >> 2)*(levelcellsperdim) + ((pz+zz) >> 2)*(levelcellsperdim2);
                    if(render[oLevel1] != 0)
                        if((d.x-xx*levelroot.w)*(d.x-xx*levelroot.w) +
                                (d.y-yy*levelroot.w)*(d.y-yy*levelroot.w) + 
                                (d.z-zz*levelroot.w)*(d.z-zz*levelroot.w) < d.w
                        )

                        {
                            unsigned int position = atomic_add(&render[oLevel1 + MAX_GPU_MAINCELLCOUNT], 1);
                            int o = ((px+xx) & 0x3) + (((py+yy) & 0x3)<< 2) + (((pz+zz) & 0x3) << 4);
                            if(bWrite)
                                stack[position + localConstants.stackOffset] = atomId + (o << 24);
                        }
                }
            }
        }
    }
}


void RenderTopLevelAtoms(
        local int* render,
        const float4 level2root,
        const int level2cellsperdimmask,
        const int level2cellsperdim,
        const int threadIdx,
        const pGPUConstants localConstants,
        global const float4* atoms
)
{
    for(int atomId = threadIdx ; atomId < localConstants.atomCnt; atomId += GPU_BLOCK_SIZE)
    {
        float4 a = atoms[atomId + localConstants.atomsOffset];
        float4 d;
        d.xyz = a.xyz - level2root.xyz;
        d.w = a.w + level2root.w*SQRT3*0.5f;
        int px = d.x/level2root.w;
        int py = d.y/level2root.w;
        int pz = d.z/level2root.w;
        d.x -= px*level2root.w;
        d.y -= py*level2root.w;
        d.z -= pz*level2root.w;
        d.w *=d.w;
        // assume grid size is larger than radius
        int xxp = (px < level2cellsperdimmask && level2root.w < d.x + a.w);
        int xxm = -(px > 0 && a.w > d.x);
        int yyp = (py < level2cellsperdimmask && level2root.w < d.y + a.w);
        int yym = -(py > 0 && a.w > d.y);
        int zzp = (pz < level2cellsperdimmask && level2root.w < d.z + a.w);
        int zz = -(pz > 0 && a.w > d.z);
        d.xyz -= 0.5f*level2root.w;
        for(; zz <= zzp; zz++)
        {
            for(int yy = yym; yy <= yyp; yy++)
            {
                for(int xx = xxm; xx <= xxp; xx++)
                {
                    if((d.x-xx*level2root.w)*(d.x-xx*level2root.w) +
                            (d.y-yy*level2root.w)*(d.y-yy*level2root.w) + 
                            (d.z-zz*level2root.w)*(d.z-zz*level2root.w) < d.w
                    )
                    {
                        int oLevel1 = (px+xx) + (py+yy)*level2cellsperdim;
                        atomic_or(&render[oLevel1], 1 << (pz+zz)) ;
                    }
                }
            }
        }
    }
}


void ComputeTopLevelWaterIntegrals(
        global float* tempIntegrals,
        local int* level3scratch,
        local float* scratch,
        const float4 level1root,
        const float threshold,
        const float level1waterW3,
        const float level1precisewaterW3,
        const int level1cellsperdimmask,
        const int threadIdx,
        const pGPUConstants localConstants,
        const int integrationBlockId,
        const int integrationBlockCount,
        global const float4* atoms,
        local const int* render
)
{
    float3 localCell = level1root.xyz;
    int threadwarpid = threadIdx & 0x1f;
    localCell.z += 0.5f * level1root.w;
    localCell.x += (0.5f + (float)(threadwarpid & level1cellsperdimmask)) * level1root.w;
    threadwarpid >>= localConstants.maincelllevel;
    localCell.y += (0.5f + (float)(threadwarpid & level1cellsperdimmask)) * level1root.w;
    if(localConstants.maincelllevel == 2)
    {
        localCell.z += ((threadwarpid >> 2) & level1cellsperdimmask) * level1root.w;
    }
    float3 preciseCell = level1root.xyz;
    preciseCell.x += (0.5f*0.125f + 0.125f*(threadIdx & 0x7)) * level1root.w;
    preciseCell.y += (0.5f*0.125f + 0.125f*((threadIdx & 0x1f) >> 3)) * level1root.w;
    preciseCell.z += (0.5f*0.125f) * level1root.w;
    for(int atomId =(threadIdx / WARP_SIZE)+integrationBlockId*GPU_BLOCK_SIZE/WARP_SIZE; atomId <localConstants.atomCnt; atomId += integrationBlockCount*GPU_BLOCK_SIZE/WARP_SIZE) //switch to atomic queue - save some mus.
    {
        threadwarpid = threadIdx & 0x1f;
        int scratchoffset = ((threadIdx/32) << 4); // 16 registers to encode up to 512 bits. each bit is set if corresponding water cell needs precise integral
        if(threadwarpid < 16)
            level3scratch[scratchoffset + threadwarpid] = 0;

        mem_fence(CLK_LOCAL_MEM_FENCE);

        float4 a = atoms[atomId + localConstants.atomsOffset];
        float sum = 0;
        float3 dd = localCell.xyz - a.xyz;
        float d2x = dd.x*dd.x;
        if(localConstants.maincelllevel == 2)
        {
            if(render[threadwarpid] == 0 )
            {
                float d = d2x + dd.y*dd.y + dd.z*dd.z;
                if(d >= threshold )
                    sum += level1waterW3/ (d*d*d);
                else
                {
                    atomic_or(&level3scratch[scratchoffset], 1 << threadwarpid);
                }
            }

            if( render[threadwarpid + WARP_SIZE] == 0)
            {
                dd.z += 0.5f*localConstants.root.w;
                float d = d2x + dd.y*dd.y + dd.z*dd.z;
                if(d >= threshold )
                    sum += level1waterW3/ (d*d*d);
                else
                {
                    atomic_or(&level3scratch[scratchoffset + 1], 1 << threadwarpid);
                }
            }
        } else
        {
            // 3
            float d1y = d2x + dd.y*dd.y;
            float d2y = d2x + (dd.y + 0.5f*localConstants.root.w)*(dd.y + 0.5f*localConstants.root.w);
            for(; threadwarpid < 512; threadwarpid+=64, dd.z += level1root.w)
            {
                if(render[threadwarpid] == 0)
                {
                    float d = d1y + dd.z*dd.z;

                    if(d >= threshold )
                        sum += level1waterW3/ (d*d*d);
                    else
                    {
                        atomic_or(&level3scratch[scratchoffset + (threadwarpid >> 5)], 1 << (threadwarpid & 31));
                    }
                }
                if( render[threadwarpid + WARP_SIZE] == 0)
                {
                    float d = d2y + dd.z*dd.z;
                    if(d >= threshold )
                        sum += level1waterW3/ (d*d*d);
                    else
                    {
                        atomic_or(&level3scratch[scratchoffset + (threadwarpid >> 5) + 1], 1 << (threadwarpid & 31));
                    }
                }
            }
        }
        mem_fence(CLK_LOCAL_MEM_FENCE);

        // have all apropriate bits set. do precise integration over all set bits for all maincells
        int precisecellid = 0;
        uint mask = level3scratch[scratchoffset];
        while (precisecellid < localConstants.maincellcount)
        {
            if(mask == 0)
            {
                precisecellid += 32;
                mask = level3scratch[scratchoffset + (precisecellid >> 5)];
                continue;
            }
            uint setbit = 31 - clz(mask);
            mask -= 1 << setbit;
            setbit += precisecellid;
            float3 dd = preciseCell.xyz - a.xyz;
            dd.x += (setbit & level1cellsperdimmask) * level1root.w;
            setbit >>= localConstants.maincelllevel;
            dd.y += (setbit & level1cellsperdimmask) * level1root.w;
            setbit >>= localConstants.maincelllevel;
            dd.z += (setbit & level1cellsperdimmask) * level1root.w;
            float dy1 = dd.y*dd.y + dd.x*dd.x;
            float dy2;
            if(WARP_SIZE == 32)
            {
                dy2 = (dd.y + 0.5f*level1root.w)*(dd.y + 0.5f*level1root.w) + dd.x*dd.x;
            }
            for(int zz = 0; zz < 8; zz++, dd.z += 0.125f*level1root.w)
            {
                float d = dy1 + dd.z*dd.z;
                sum += level1precisewaterW3 / (d*d*d);
                if(WARP_SIZE == 32)
                {
                    d = dy2 + dd.z*dd.z;
                    sum += level1precisewaterW3 / (d*d*d);
                }
            }

        }

        // all integrals for an atom completed, do a stable reduction for sum.
        int scratchthread = threadIdx >> 2;
        int scratchmask = threadIdx & 0x3;
        if(scratchmask == 3)
        {
            scratch[scratchthread] = sum;
        }
        mem_fence(CLK_LOCAL_MEM_FENCE);
        if(scratchmask == 1)
        {
            scratch[scratchthread] += sum;
        }
        mem_fence(CLK_LOCAL_MEM_FENCE);
        if(scratchmask == 2)
        {
            scratch[scratchthread] += sum;
        }
        mem_fence(CLK_LOCAL_MEM_FENCE);
        if(scratchmask == 0)
        {
            scratchmask = threadIdx & (WARP_SIZE - 1);// we have threads 0,4,8,12,16,20,24,28 left
            scratch[scratchthread] += sum;
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if(scratchmask < 16)
            {
                scratch[scratchthread] += scratch[scratchthread + 4];
            }
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if(scratchmask < 8)
            {
                scratch[scratchthread] += scratch[scratchthread + 2];
            }
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if(scratchmask == 0)
            {
                tempIntegrals[atomId + localConstants.atomsOffset*integrationBlockCount + integrationBlockId*localConstants.atomPitch] += scratch[scratchthread] + scratch[scratchthread + 1];
            }
        }

    }
}


void RenderCellAtomsAndWater(
        local int* render,
        local float4* level3,
        local int* level3scratch,
        local int* level1AtomIdId,
        const int level1WaterIdId,
        const int level1AtomIdCnt,
        const int level1WaterIdCnt,
        const float4 level1root,
        const float lastLevelW,
        const float lastLevelW3,
        const int threadIdx,
        const pGPUConstants localConstants,
        global const int* stack,
        global const float4* atoms,
        global const float4* water,
        local const int* ss,
        local const int* vv
)
{
    float px = level1root.x + ( (float)(threadIdx&(LAST_LEVEL_GRID_SIZE-1)))*lastLevelW;

    while (*level1AtomIdId <level1AtomIdCnt) // render atoms
    {
        if(threadIdx % WARP_SIZE == 0)
        {
            //load next atom
            int position = atomic_add(level1AtomIdId, 1);
            if(position >=level1AtomIdCnt)
            {
                level3scratch[threadIdx / WARP_SIZE ] = -1;
            }
            else
            {
                level3scratch[threadIdx / WARP_SIZE ] = stack[position + localConstants.stackOffset];
            }
        }
        mem_fence(CLK_LOCAL_MEM_FENCE);
        if(level3scratch[threadIdx / WARP_SIZE ] < 0) break;

        float4 a = atoms[level3scratch[threadIdx / WARP_SIZE ] + localConstants.atomsOffset];

        float dx = px - a.x;
        dx=a.w*a.w - dx*dx;
        if(dx > 0)
        {
            int py = (a.y-a.w - level1root.y)/lastLevelW;
            if(py < 0) py = 0;
            int py2 = (a.y+a.w - level1root.y)/lastLevelW;
            if(py2 > LAST_LEVEL_GRID_SIZE -1) py2 = LAST_LEVEL_GRID_SIZE -1;
            float dy = level1root.y +py*lastLevelW - a.y;
            for(; py <= py2; py++, dy+= lastLevelW)
            {
                float dxy = dx - dy*dy;
                if(dxy <= 0) continue;
                float dz = sqrt(dxy);
                int dz1 = (int)((a.z - level1root.z - dz)/lastLevelW);
                int dz2 = (int)((a.z - level1root.z + dz)/lastLevelW);
                if(dz1 >= LAST_LEVEL_GRID_SIZE - 1) continue;
                if(dz2 <= 0) continue;
                unsigned int mask = 0;
                if(dz2 < LAST_LEVEL_GRID_SIZE)
                {
                    mask = 0xffffffff - ((2 << dz2) - 1);
                }
                if(dz1 >= 0)
                {
                    mask |= ((2 << dz1) - 1);
                }
                atomic_and(&render[py*LAST_LEVEL_GRID_SIZE + (threadIdx&(LAST_LEVEL_GRID_SIZE-1))], mask);
            }
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // render water
    for(int wPosId = level1WaterIdId+ threadIdx / WARP_SIZE2; wPosId < level1WaterIdCnt; wPosId += GPU_BLOCK_SIZE_KERNEL2/WARP_SIZE2)
    {
        int aId = stack[wPosId  + localConstants.stackOffset];
        float4 a = water[aId + localConstants.waterOffset];

        int px = (threadIdx & 0xf) - 8 + (a.x-level1root.x)/lastLevelW;
        float dx = a.x - level1root.x - px*lastLevelW;
        dx=a.w*a.w - dx*dx;
        if(dx > 0 && px >=0 && px <LAST_LEVEL_GRID_SIZE)
        {
            int py = (a.y-a.w - level1root.y)/lastLevelW;
            int py2 = (a.y+a.w - level1root.y)/lastLevelW;
            if(py < 0) py = 0;
            if(py2 > LAST_LEVEL_GRID_SIZE -1) py2 = LAST_LEVEL_GRID_SIZE -1;
            float dy = level1root.y +py*lastLevelW - a.y;
            for(; py <= py2; py++, dy+= lastLevelW) //using hardcoded 16x1 grid, water size can't be greater than 1.75f! 1 warp processes 2 or 4 waters at once.
            {
                float dxy = dx - dy*dy;
                if(dxy <= 0) continue;
                float dz = sqrt(dxy);
                int dz1 = (int)((a.z - level1root.z - dz)/lastLevelW);
                int dz2 = (int)((a.z - level1root.z + dz)/lastLevelW);
                if(dz1 >= LAST_LEVEL_GRID_SIZE - 1) continue;
                if(dz2 <= 0) continue;
                unsigned int mask = 0xffffffff;
                if(dz2 < LAST_LEVEL_GRID_SIZE)
                {
                    mask = (2 << dz2) - 1;
                }
                if(dz1 >= 0)
                {
                    mask &= 0xffffffff - ((2 << dz1) - 1);
                }
                atomic_or(&render[py*LAST_LEVEL_GRID_SIZE + px], mask);
            }
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);


    if(threadIdx < LEVEL3_SIZE*2)
    {
        level3scratch[threadIdx] = 0;
        level3scratch[threadIdx + LEVEL3_SIZE*2] = 0;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // low level render complete. generate top level data level3 is 4*4*4. each 8 threads compute 1 cell. each thread has fixed x and loops over all y, all Z are processed in bulk
    // this takes 30mus - 15% of render time. atomic_adds may be optimized for a small ammount of 8mus
    for(int id = threadIdx; id < 512; id+=GPU_BLOCK_SIZE_KERNEL2)
    {
        int x = id & 0x1f;
        int z = (id>>2) & 0x18;
        int y = (id>>4) & 0x18;
        int vy = 0;
        int vz = 0;
        int vs = 0;
        int maskid = y*LAST_LEVEL_GRID_SIZE + x;
        for(int yy = 0; yy < 8; yy++, maskid+=LAST_LEVEL_GRID_SIZE)
        {
            unsigned int mask = render[maskid];
            mask >>= z;
            int s = ss[mask & 0xf];
            vs+=s;
            vy+=s*yy;
            vz+=vv[mask & 0xf];
            mask = (mask >> 4) & 0xf;
            s = ss[mask];
            vs+=s;
            vy+=s*yy;
            vz+=vv[mask] + 4*s;
        }
        // level3
        int offset3 = (x >> 3) + (y >> 1) + (z << 1);

        atomic_add(&(level3scratch[offset3 + 0*LEVEL3_SIZE]), vs); // use in-warp reduction - save 8mus. no atomic_add for vector or vector elements!
        atomic_add(&(level3scratch[offset3 + 1*LEVEL3_SIZE]), vs*x);
        atomic_add(&(level3scratch[offset3 + 2*LEVEL3_SIZE]), vy + vs*y);
        atomic_add(&(level3scratch[offset3 + 3*LEVEL3_SIZE]), vz + vs*z);
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if(threadIdx < LEVEL3_SIZE) //10mus for this
    {
        level3[threadIdx].w = level3scratch[threadIdx + 0*LEVEL3_SIZE]*lastLevelW3;
        level3[threadIdx].x = lastLevelW*(float)level3scratch[threadIdx + 1*LEVEL3_SIZE]/(float)level3scratch[threadIdx + 0*LEVEL3_SIZE] + level1root.x;
        level3[threadIdx].y = lastLevelW*(float)level3scratch[threadIdx + 2*LEVEL3_SIZE]/(float)level3scratch[threadIdx + 0*LEVEL3_SIZE] + level1root.y;
        level3[threadIdx].z = lastLevelW*(float)level3scratch[threadIdx + 3*LEVEL3_SIZE]/(float)level3scratch[threadIdx + 0*LEVEL3_SIZE] + level1root.z;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
}


void ComputeCoarseIntegral(
        local int* level1IntegIdCnt,
        local int* level3scratch,
        global int* atomIdBuffer,
        global float* tempIntegrals,
        const int threadIdx,
        const float threshold,
        const int atomIdBufferIdx,
        const int integrationBlockOffset,
        const pGPUConstants localConstants,
        global const float4* atoms,
        local const float4* level3
)
{
    for(int atomId = threadIdx ; atomId < localConstants.atomCnt; atomId += GPU_BLOCK_SIZE_KERNEL2)
    {
        float sum = 0;
        float4 a = atoms[atomId + localConstants.atomsOffset];
        int cellId = 0;
        int mask =0;
        while( cellId < 32)
        {
            if(level3[cellId].w > 0)
            {
                float3 v = a.xyz - level3[cellId].xyz;
                float d = v.x*v.x + v.y*v.y + v.z*v.z;
                if(d >= threshold )
                    sum += level3[cellId].w/ (d*d*d);

                else
                {
                    mask |= 1 << cellId;    
                }

            }
            cellId ++;
        }
        int mask2 =0;
        while( cellId < 64)
        {
            if(level3[cellId].w > 0)
            {
                float3 v = a.xyz - level3[cellId].xyz;
                float d = v.x*v.x + v.y*v.y + v.z*v.z;
                if(d >= threshold )
                    sum += level3[cellId].w/ (d*d*d);

                else
                {
                    mask2 |= 1 << (cellId - 32);    
                }

            }
            cellId ++;
        }
        tempIntegrals[atomId + integrationBlockOffset] +=sum;                
        if(mask || mask2)
        {
            int position = 4*atomic_add(level1IntegIdCnt, 1);
            if(position < LEVEL3_SIZE*4)
            {
                level3scratch[position] = atomId;
                level3scratch[position+1] = mask;
                level3scratch[position+2] = mask2;
            }
            else
            {
                atomIdBuffer[atomIdBufferIdx + position] = atomId;
                atomIdBuffer[atomIdBufferIdx + position + 1] = mask;
                atomIdBuffer[atomIdBufferIdx + position + 2] = mask2;
            }
        }
    }
}

/*
void ComputePreciseIntegral(
    global float* tempIntegrals,    
    local float* scratch,
    local int* level3scratch,
    global const int* atomIdBuffer,
    const int level1IntegIdCnt,
    const int threadIdx, 
    const int atomIdBufferIdx,
    const int integrationBlockOffset,
    const pGPUConstants localConstants, 
    const float4 level1root,
    const float lastLevelW, 
    const float lastLevelW3, 
    local const int* render, 
    global const float4* atoms
)
{
}
 */

__kernel void
PowerBornKernel1(__global const float4* atoms, 
        __global float4* water,
        __global pGPUConstants* constants,
        __global float* result,
        __global float* resultArea,
        __global const float4* waterSeeds,
        __global const int* waterMask,
        __global int* stack,
        const int integrationBlockCount
)
{
    __local int ss[16];
    __local int vv[16];
    __local int render[LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE];
    __local int level3scratch[LEVEL3_SIZE*4 + 1];
    __local int stackPos; // 

    __local unsigned int level1cellsperdim ;
    __local unsigned int level1cellsperdim2 ;
    __local unsigned int level1cellsperdimmask ;
    __local float4 level1root;
    __local float cellOffset  ;


    __local unsigned int level2cellsperdim ;
    __local unsigned int level2cellsperdim2 ;
    __local unsigned int level2cellsperdimmask ;
    __local float4 level2root;

    __local pGPUConstants localConstants;

    __local float threshold;
    __local int level1CellDataOffset;
    __local int level1CellDataTotal;
    __local int level1WaterIdCnt;
    __local int level1IntegIdCnt;
    __local float level1waterW3;
    __local float level1precisewaterW3;
    __local int level1CellCnt;
    __local int level2PreciseWaterCnt;
    __local float level2waterW3;
    __local int largelevel2cellsperdimmask;
    __local int curStackPos;
    int threadIdx = get_local_id(0); // this is 20mus faster than calling get_local_id(0) all the time!
    // initial population of coarse grid
    {
        // read init data
        if(threadIdx == 0)
        {
            int blockIdx = get_group_id(0);
            localConstants = constants[blockIdx];

            level1cellsperdim = 1 << localConstants.maincelllevel;
            level1cellsperdimmask = level1cellsperdim - 1;
            cellOffset  =  - ((float)level1cellsperdim / 2);
            threshold = 0.5f*localConstants.root.w / level1cellsperdim;
            threshold*= threshold*localConstants.integration_factor;
            level1cellsperdim2 = level1cellsperdim *level1cellsperdim;
            level1root.w = localConstants.root.w / level1cellsperdim;
            level1root.xyz = localConstants.root.xyz + (cellOffset)*level1root.www;
            cellOffset += 0.5f;
            level1waterW3 = level1root.w*level1root.w*level1root.w;
            level1precisewaterW3 = level1waterW3 / 512.0f;

            level2cellsperdim = level1cellsperdim;
            level2cellsperdimmask = level1cellsperdimmask;
            level2cellsperdim2 = level1cellsperdim2;
            level2waterW3 = level1waterW3;
            level2root = level1root;
            if(localConstants.maxLevel > 8)
            {
                level2cellsperdim <<=  2;
                level2cellsperdimmask = level2cellsperdim -1;
                level2cellsperdim2 <<= 4;
                level2waterW3 = level2waterW3/ 64.0f;
                level2root.w *= 0.25f;
            }
            level1WaterIdCnt = 0;
        }
        else if(threadIdx >= WARP_SIZEL)
        {
            for(int atomid = threadIdx - WARP_SIZEL; atomid < LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE; atomid += GPU_BLOCK_SIZE - WARP_SIZEL) 
                render[atomid] = 0;
        }

        barrier(CLK_LOCAL_MEM_FENCE); 

        for(int atomid = threadIdx; atomid < localConstants.atomPitch*integrationBlockCount; atomid += GPU_BLOCK_SIZE)
        {
            result[atomid + localConstants.atomsOffset*integrationBlockCount] = 0;
        }

        const unsigned int blockIdx = get_group_id(0);

        for(int atomid = threadIdx; atomid < localConstants.atomCnt; atomid += GPU_BLOCK_SIZE)
        {
            float4 a = atoms[atomid + localConstants.atomsOffset];
            int maskadr = (atomid + localConstants.atomsOffset)*(WATERMASK_SIZE/32);
            int wid = 0;
            unsigned int mask = waterMask[maskadr];
            int waterSeedsCnt; 
            int waterSeedsOffset;
            float area = 0; 
            if(a.w < WATERMASK_THRESHOLD)
            {
                waterSeedsCnt = WATERMASK_SMALL;
                waterSeedsOffset = 0;
            }
            else
            {
                waterSeedsCnt = WATERMASK_LARGE;
                waterSeedsOffset = WATERMASK_SMALL;
            }
            while (wid < waterSeedsCnt)
            {
                if(mask == 0)
                {
                    wid += 32;
                    maskadr +=1;
                    mask = waterMask[maskadr];
                    continue;
                }
                uint setbit = 31 - clz(mask);
                mask -= 1 << setbit;
                setbit += wid;
                if(setbit < waterSeedsCnt)
                {
                    float4 w = waterSeeds[setbit + waterSeedsOffset];
                    area += w.w;
                    w.xyz = a.xyz + a.www*w.xyz;
                    w.w = 1.4f;
                    int pos = atomic_add(&level1WaterIdCnt,1);
                    water[pos + localConstants.waterOffset] = w;
                }
            }
            resultArea[atomid + localConstants.atomsOffset] = area * a.w * a.w;
        } 

        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE); 

        // count atoms and water in top level cells
        // out stack - not modified since bWrite == 0, atoms and waters are only counted here
        // in/out render - render will contain counts for atoms and water atoms per top level cell
        if(localConstants.maxLevel > 8)
        {
            AssignTopLevelAtomAndWaterToCellsLargeGrid(
                    stack, render, 0,
                    level2root, level2cellsperdimmask, level1cellsperdim, level1cellsperdim2, threadIdx, localConstants, atoms, water
            );
        }
        else
        {
            AssignTopLevelAtomAndWaterToCells(
                    stack, render, 0,
                    level1root, level1cellsperdimmask, level1cellsperdim, level1cellsperdim2, threadIdx, localConstants, atoms, water
            );
        }


        barrier(CLK_LOCAL_MEM_FENCE); 
        // counts computed, allocate memory - do a reduction to compute running sum. Reusing buffers and supporting multiple resolutions, so can't easily write a loop here.
        for(int id = threadIdx; id < localConstants.maincellcount; id +=GPU_BLOCK_SIZE)
        {
            if(render[id] > 0)
                render[id + MAX_GPU_MAINCELLCOUNT] += render[id] + 4;
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if((id & 0x1) == 0)
                render[id + MAX_GPU_MAINCELLCOUNT] += render[id + MAX_GPU_MAINCELLCOUNT + 1];
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if((id & 0x3) == 0)
                level3scratch[(id >>2) + 128] = render[id + MAX_GPU_MAINCELLCOUNT] + render[id + MAX_GPU_MAINCELLCOUNT + 2];
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if((id & 0x7) == 0)
                level3scratch[(id >>3) + 64] = level3scratch[(id >>2) + 128] + level3scratch[(id >>2) + 129];
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if((id & 0xf) == 0)
                level3scratch[(id >>4) + 32] = level3scratch[(id >>3) + 64] + level3scratch[(id >>3) + 65];
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if((id & 0x1f) == 0)
                level3scratch[(id >>5) + 16] = level3scratch[(id >>4) + 32] + level3scratch[(id >>4) + 33];
            mem_fence(CLK_LOCAL_MEM_FENCE);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if(localConstants.maincellcount == 64)
        {
            if(threadIdx == 0)
            {
                level3scratch[8] = level3scratch[16] + level3scratch[17];
                stackPos = level3scratch[8] + MAX_GPU_MAINCELLCOUNT;
                //                octTree[0+ octTreeOffset].celltype[0] = stackPos;
            }
        } 
        else if(threadIdx < 8)
        {
            level3scratch[threadIdx + 8] = level3scratch[threadIdx*2 + 16] + level3scratch[threadIdx*2 + 17];
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if(threadIdx < 4)
                level3scratch[threadIdx + 4] = level3scratch[threadIdx*2 + 8] + level3scratch[threadIdx*2 + 9];
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if(threadIdx < 2)
                level3scratch[threadIdx + 2] = level3scratch[threadIdx*2 + 4] + level3scratch[threadIdx*2 + 5];
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if(threadIdx == 0)
            {
                level3scratch[1] = level3scratch[2] + level3scratch[3];
                stackPos = level3scratch[1] + MAX_GPU_MAINCELLCOUNT;
                stack[0+ localConstants.stackOffset] = stackPos;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); 

        for(int id = threadIdx; id < localConstants.maincellcount; id +=GPU_BLOCK_SIZE)
        {
            if(render[id] > 0)
            {
                //compute stack alloc 
                int maincellPos = MAX_GPU_MAINCELLCOUNT;
                if(id & 0x1)
                {
                    maincellPos += render[id + MAX_GPU_MAINCELLCOUNT - 1] - render[id + MAX_GPU_MAINCELLCOUNT];
                }
                if(id & 0x2)
                {
                    maincellPos += render[(id & 0xfffffc) + MAX_GPU_MAINCELLCOUNT];
                }
                int mainstep = 128;
                int binId = id >> 2;
                while (mainstep > 0 && binId > 0)
                {
                    if(binId & 1)
                        maincellPos += level3scratch[(binId & 0xffe) + mainstep];
                    binId >>= 1;
                    mainstep >>= 1;
                };
                //maincellPos contains offset of cell data from the start of stack - all is ordered by cellid
                int atomstart = maincellPos + 4;
                int waterstart = atomstart + render[id];
                int waterend = maincellPos + render[id + MAX_GPU_MAINCELLCOUNT];
                if((id & 0x1) == 0)
                {
                    waterend -= render[id + MAX_GPU_MAINCELLCOUNT + 1];
                }
                stack[maincellPos + localConstants.stackOffset] = atomstart;
                stack[maincellPos + 1 + localConstants.stackOffset] = waterstart;
                stack[maincellPos + 2 + localConstants.stackOffset] = waterend;
                stack[maincellPos + 3 + localConstants.stackOffset] = id;
                //                octTree[threadIdx + octTreeOffset].celltype[1] = atomstart;
                //                octTree[threadIdx + octTreeOffset].celltype[2] = waterstart;
                render[id + MAX_GPU_MAINCELLCOUNT] = waterstart;
                render[id] = atomstart;
            }
            mem_fence(CLK_LOCAL_MEM_FENCE);
            if(render[id] == 0)
            {
                render[id + MAX_GPU_MAINCELLCOUNT] = 0;

            }
            mem_fence(CLK_LOCAL_MEM_FENCE);
        };

        barrier(CLK_LOCAL_MEM_FENCE); 
        // assign atoms to cells
        // should produce exactly the same result as the code that computes counts
        // fills stack with atom and water ids that belong to particular cell
        // out stack
        // in/out render - render contains offsets this function needs to write atomids to. each write increments the offset and at the end the array contains final locations that were written
        // one can compare these locations with what stack[] has at maincell+0 maincell+1 for this and next cell for debug purposes 
        if(localConstants.maxLevel > 8)
        {
            AssignTopLevelAtomAndWaterToCellsLargeGrid(
                    stack, render, 1,
                    level2root, level2cellsperdimmask, level1cellsperdim, level1cellsperdim2, threadIdx, localConstants, atoms, water
            );
        }
        else
        {
            AssignTopLevelAtomAndWaterToCells(
                    stack, render, 1,
                    level1root, level1cellsperdimmask, level1cellsperdim, level1cellsperdim2, threadIdx, localConstants, atoms, water
            );
        }

        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE); 

        // save mid level render data
        for(int atomid = threadIdx; atomid < MAX_GPU_MAINCELLCOUNT; atomid += GPU_BLOCK_SIZE ) 
            stack[atomid + localConstants.stackOffset - 3*MAX_GPU_MAINCELLCOUNT] = render[atomid];
    }

    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE); 

    if(localConstants.maxLevel > 8) //large grid - sort atoms on lower levels
    {
        if(threadIdx == 0)
        {
            level2PreciseWaterCnt = 0;
            level1CellCnt = 0;
            level1CellDataTotal = stackPos;
            level1CellDataOffset = MAX_GPU_MAINCELLCOUNT;
            largelevel2cellsperdimmask = level2cellsperdimmask & 0xffc;
            level1IntegIdCnt = 0;
            threshold /= 16.0f;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        while (level1CellDataOffset < level1CellDataTotal)
        {
            // read index
            barrier(CLK_LOCAL_MEM_FENCE);
            if(threadIdx == 0)
            {
                level1IntegIdCnt = 0;
                while (level1IntegIdCnt < MAX_GPU_MAINCELLCOUNT/64 && level1CellDataOffset < level1CellDataTotal)
                {
                    int atomStart = stack[level1CellDataOffset + localConstants.stackOffset]; // can do it as 1 16-byte read.
                    int waterStart = stack[level1CellDataOffset + 1 + localConstants.stackOffset];
                    int waterEnd = stack[level1CellDataOffset + 2 + localConstants.stackOffset];
                    int cellId = stack[level1CellDataOffset + 3 + localConstants.stackOffset];
                    ss[level1IntegIdCnt] = atomStart;
                    vv[level1IntegIdCnt] = waterStart;
                    vv[level1IntegIdCnt + MAX_GPU_MAINCELLCOUNT/64] = cellId;
                    level1CellDataOffset = waterEnd;
                    ss[level1IntegIdCnt + 1] = waterEnd + 4;
                    level1IntegIdCnt++;
                }
                level1CellCnt += level1IntegIdCnt;
                curStackPos = stackPos;
            }
            for(int atomid = threadIdx ; atomid < LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE; atomid += GPU_BLOCK_SIZE )
            {
                render[atomid] = 0;
            }
            barrier(CLK_LOCAL_MEM_FENCE);

            if(threadIdx < level1IntegIdCnt)
            {
                int cellId = vv[threadIdx + MAX_GPU_MAINCELLCOUNT/64];
                stack[level1CellCnt - level1IntegIdCnt + threadIdx + localConstants.stackOffset] = cellId;
                vv[threadIdx + MAX_GPU_MAINCELLCOUNT/64] =(
                        ((cellId << 6) & (largelevel2cellsperdimmask << (2*localConstants.maincelllevel + 4))) +
                        ((cellId << 4) & (largelevel2cellsperdimmask << (localConstants.maincelllevel + 2))) +
                        ((cellId << 2) & (largelevel2cellsperdimmask)) 
                );
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            int slot = 0;
            for(int atomIdId = ss[0] + threadIdx; atomIdId < level1CellDataOffset; atomIdId+=GPU_BLOCK_SIZE)
            {
                int atomId = stack[atomIdId + localConstants.stackOffset];
                while (atomIdId >= ss[slot+1]) slot++;
                if(atomIdId >= ss[slot+1] - 4) continue; // in index space
                int cell = slot*64 + (atomId >> 24);
                if(atomIdId >= vv[slot]) cell += MAX_GPU_MAINCELLCOUNT; // water atom
                atomic_add(&render[cell], 1);
            }
            barrier(CLK_LOCAL_MEM_FENCE);

            for(int id = threadIdx; id < MAX_GPU_MAINCELLCOUNT; id +=GPU_BLOCK_SIZE)
            {
                if(render[id] > 0)
                    render[id + MAX_GPU_MAINCELLCOUNT] += render[id] + 4;
                else
                    render[id + MAX_GPU_MAINCELLCOUNT] = 0;
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if((id & 0x1) == 0)
                    render[id + MAX_GPU_MAINCELLCOUNT] += render[id + MAX_GPU_MAINCELLCOUNT + 1];
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if((id & 0x3) == 0)
                    level3scratch[(id >>2) + 128] = render[id + MAX_GPU_MAINCELLCOUNT] + render[id + MAX_GPU_MAINCELLCOUNT + 2];
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if((id & 0x7) == 0)
                    level3scratch[(id >>3) + 64] = level3scratch[(id >>2) + 128] + level3scratch[(id >>2) + 129];
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if((id & 0xf) == 0)
                    level3scratch[(id >>4) + 32] = level3scratch[(id >>3) + 64] + level3scratch[(id >>3) + 65];
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if((id & 0x1f) == 0)
                    level3scratch[(id >>5) + 16] = level3scratch[(id >>4) + 32] + level3scratch[(id >>4) + 33];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if(threadIdx < 8)
            {
                level3scratch[threadIdx + 8] = level3scratch[threadIdx*2 + 16] + level3scratch[threadIdx*2 + 17];
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if(threadIdx < 4)
                    level3scratch[threadIdx + 4] = level3scratch[threadIdx*2 + 8] + level3scratch[threadIdx*2 + 9];
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if(threadIdx < 2)
                    level3scratch[threadIdx + 2] = level3scratch[threadIdx*2 + 4] + level3scratch[threadIdx*2 + 5];
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if(threadIdx == 0)
                {
                    level3scratch[1] = level3scratch[2] + level3scratch[3];
                    stackPos += level3scratch[1];
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE); 

            for(int id = threadIdx; id < 64*level1IntegIdCnt; id +=GPU_BLOCK_SIZE)
            {
                if(render[id] > 0)
                {
                    int maincellPos = curStackPos;
                    //compute stack alloc 
                    if(id & 0x1)
                    {
                        maincellPos += render[id + MAX_GPU_MAINCELLCOUNT - 1] - render[id + MAX_GPU_MAINCELLCOUNT];
                    }
                    if(id & 0x2)
                    {
                        maincellPos += render[(id & 0xfffffc) + MAX_GPU_MAINCELLCOUNT];
                    }
                    int mainstep = 128;
                    int binId = id >> 2;
                    while (mainstep > 0 && binId > 0)
                    {
                        if(binId & 1)
                            maincellPos += level3scratch[(binId & 0xffe) + mainstep];
                        binId >>= 1;
                        mainstep >>= 1;
                    };
                    //maincellPos contains offset of cell data from the start of stack - all is ordered by cellid
                    int atomstart = maincellPos + 4;
                    int waterstart = atomstart + render[id];
                    int waterend = maincellPos + render[id + MAX_GPU_MAINCELLCOUNT];
                    if((id & 0x1) == 0)
                    {
                        waterend -= render[id + MAX_GPU_MAINCELLCOUNT + 1];
                    }
                    stack[maincellPos + localConstants.stackOffset] = atomstart;
                    stack[maincellPos + 1 + localConstants.stackOffset] = waterstart;
                    stack[maincellPos + 2 + localConstants.stackOffset] = waterend;
                    stack[maincellPos + 3 + localConstants.stackOffset] = vv[(id >> 6) + MAX_GPU_MAINCELLCOUNT/64] + (id & 0x3) + ((id & 0xc) << localConstants.maincelllevel) + ((id & 0x30) << (2*localConstants.maincelllevel)); // need correct cellid here
                    render[id + MAX_GPU_MAINCELLCOUNT] = waterstart;
                    render[id] = atomstart;
                }
                mem_fence(CLK_LOCAL_MEM_FENCE);
                if(render[id] == 0)
                {
                    render[id + MAX_GPU_MAINCELLCOUNT] = 0;
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            slot = 0;
            for(int atomIdId = ss[0] + threadIdx; atomIdId < level1CellDataOffset; atomIdId+=GPU_BLOCK_SIZE)
            {
                int atomId = stack[atomIdId + localConstants.stackOffset];
                while (atomIdId >= ss[slot+1]) slot++;
                if(atomIdId >= ss[slot+1] - 4) continue; // in index space
                int cell = slot*64 + (atomId >> 24);
                if(atomIdId >= vv[slot]) cell += MAX_GPU_MAINCELLCOUNT; // water atom
                if(render[cell] > 0)
                {
                    int pos = atomic_add(&render[cell], 1);
                    stack[pos + localConstants.stackOffset] = atomId & 0xffffff;
                }
            }

            barrier(CLK_LOCAL_MEM_FENCE);
        }


        for(int atomid = threadIdx ; atomid < LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE; atomid += GPU_BLOCK_SIZE ) 
            render[atomid] = 0;
        barrier(CLK_LOCAL_MEM_FENCE); 

        RenderTopLevelAtoms(
                render,
                level2root, level2cellsperdimmask, level2cellsperdim, threadIdx, localConstants, atoms
        );

        barrier(CLK_LOCAL_MEM_FENCE); 
        // save mid level render data
        for(int atomid = threadIdx; atomid < LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE; atomid += GPU_BLOCK_SIZE ) 
            stack[atomid + localConstants.stackOffset - LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE] = render[atomid];


        if(threadIdx == 0)
        {
            stack[localConstants.stackOffset - MAX_GPU_MAINCELLCOUNT*4 +2] = level1CellCnt;
        }
    }

    if(threadIdx == 0)
    {
        stack[localConstants.stackOffset - MAX_GPU_MAINCELLCOUNT*4] = stackPos;
        stack[localConstants.stackOffset - MAX_GPU_MAINCELLCOUNT*4 +1] = level1CellDataTotal;
    }
    barrier(CLK_LOCAL_MEM_FENCE| CLK_GLOBAL_MEM_FENCE);
}


__kernel void
PowerBornTopLevelWaterIntegral(__global const float4* atoms, 
        __global const pGPUConstants* constants,
        __global const int* stack,
        __global float* result,
        __global int* atomIdBuffer,
        const int integrationBlocksScale
)
{
    __local int render[LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE];
    __local float scratch[GPU_BLOCK_SIZE/4];
    __local int level3scratch[LEVEL3_SIZE*4 + 1];

    unsigned int level1cellsperdim ;
    unsigned int level1cellsperdim2 ;
    unsigned int level1cellsperdimmask ;
    float4 level1root;
    float cellOffset  ;


    unsigned int level2cellsperdim ;
    unsigned int level2cellsperdim2 ;
    unsigned int level2cellsperdimmask ;
    float4 level2root;

    pGPUConstants localConstants;

    float threshold;
    __local int level1WaterIdCnt;
    __local int level1IntegIdCnt;
    __local float level1waterW3;
    float level1precisewaterW3;
    __local int level1CellCnt;
    __local int level2PreciseWaterCnt;
    float level2waterW3;
    int threadIdx = get_local_id(0); // this is 20mus faster than calling get_local_id(0) all the time!
    int integrationBlockId = get_group_id(1);
    int integrationBlockCount = get_num_groups(1);
    // initial population of coarse grid
    int atomIdBufferIdx;
    int blockIdx = get_group_id(0);
    localConstants = constants[blockIdx];
    level1cellsperdim = 1 << localConstants.maincelllevel;
    level1cellsperdimmask = level1cellsperdim - 1;
    level1cellsperdim2 = level1cellsperdim *level1cellsperdim;


    atomIdBufferIdx = (get_group_id(0)*integrationBlockCount + integrationBlockId)*localConstants.atomPitch*integrationBlocksScale;
    level2cellsperdim = level1cellsperdim;
    level2cellsperdimmask = level1cellsperdimmask;
    level2cellsperdim2 = level1cellsperdim2;
    cellOffset  =  - ((float)level1cellsperdim / 2);
    threshold = 0.5f*localConstants.root.w / level1cellsperdim;
    threshold*= threshold*localConstants.integration_factor;
    level1root.w = localConstants.root.w / level1cellsperdim;
    level1root.xyz = localConstants.root.xyz + (cellOffset)*level1root.www;
    cellOffset += 0.5f;
    level1waterW3 = level1root.w*level1root.w*level1root.w;
    level1precisewaterW3 = level1waterW3 / 512.0f;
    level2waterW3 = level1waterW3;
    level2root = level1root;
    if(localConstants.maxLevel > 8)
    {
        level2cellsperdim <<=  2;
        level2cellsperdimmask = level2cellsperdim -1;
        level2cellsperdim2 <<= 4;
        level2waterW3 = level2waterW3/ 64.0f;
        level2root.w *= 0.25f;
    }

    {
        // read init data
        if(threadIdx == 0)
        {
            level1WaterIdCnt = 0;
        }

        barrier(CLK_LOCAL_MEM_FENCE); 


        // load cached top level render data
        for(int atomid = threadIdx; atomid < MAX_GPU_MAINCELLCOUNT; atomid += GPU_BLOCK_SIZE ) 
            render[atomid] = stack[atomid + localConstants.stackOffset - 3*MAX_GPU_MAINCELLCOUNT];

        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE); 

        // do coarse and precise integrals over main water cells
        // uses atom counts in render to decide if a cell is empty
        // modifies tempIntegrals - adds integral values over top level water
        // modifies level3scratch, scratch - temp data

        ComputeTopLevelWaterIntegrals(
                result, level3scratch, scratch,
                level1root, threshold, level1waterW3, level1precisewaterW3, level1cellsperdimmask, threadIdx, localConstants, integrationBlockId, integrationBlockCount, atoms, render
        );

    }

    barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE); 

    if(localConstants.maxLevel > 8) //large grid - sort atoms on lower levels
    {
        threshold /= 16.0f;
        if(threadIdx == 0)
        {
            level2PreciseWaterCnt = 0;
            level1CellCnt = stack[localConstants.stackOffset - MAX_GPU_MAINCELLCOUNT*4 + 2];
            level1IntegIdCnt = 0;
        }

        // load cached mid level render data
        for(int atomid = threadIdx; atomid < LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE; atomid += GPU_BLOCK_SIZE ) 
            render[atomid] = stack[atomid + localConstants.stackOffset - LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE];

        barrier(CLK_LOCAL_MEM_FENCE); 
        for(int atomId = threadIdx; atomId < localConstants.atomCnt; atomId += GPU_BLOCK_SIZE)
        {
            float sum = 0;
            float3 a = atoms[atomId + localConstants.atomsOffset].xyz - level1root.xyz - 0.5f*level2root.w;
            for(int cellId = integrationBlockId; cellId < level1CellCnt; cellId+=integrationBlockCount)
            {
                int cell = stack[cellId + localConstants.stackOffset];
                int3 c;
                c.x = (cell & level1cellsperdimmask) << 2;
                cell >>=localConstants.maincelllevel;
                c.y = (cell & level1cellsperdimmask) << 2;
                cell >>= localConstants.maincelllevel;
                c.z = cell << 2;
                float cpx = a.x - c.x*level2root.w;
                for(int xx= 0; xx < 4; ++xx, cpx -=level2root.w)
                {
                    float cpy = a.y - c.y*level2root.w;
                    for(int yy= 0; yy < 4; ++yy, cpy -=level2root.w)
                    {
                        float dy = cpx*cpx + cpy*cpy;
                        int adr = c.x+xx + ((c.y+yy)*level2cellsperdim);
                        int mask = render[adr] >> c.z;    
                        float cpz = a.z - c.z*level2root.w;
                        for(int zz= 0; zz < 4; ++zz, cpz -=level2root.w, mask >>= 1)
                        {
                            if((mask & 1) == 0)
                            {
                                float d = dy + cpz*cpz;
                                if(d < threshold)
                                {
                                    int pos = atomic_add(&level2PreciseWaterCnt, 1);
                                    atomIdBuffer[atomIdBufferIdx + pos] = atomId + ((adr + (c.z + zz)*level2cellsperdim2) << 16); 
                                }
                                else sum+= level2waterW3/(d*d*d);
                            }
                        }
                    }
                }
            }
            result[atomId + localConstants.atomsOffset*integrationBlockCount + localConstants.atomPitch*integrationBlockId] += sum;
        }
        // now do precise

        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
        if(level2PreciseWaterCnt > 0)
        {
            float level2rootw125 = 0.125f*level2root.w;
            float2 localCell = level2root.xy;
            localCell.x += (threadIdx & 0x7) * level2rootw125;
            localCell.y += ((threadIdx>>3) & 0x3) * level2rootw125;
            int atomIdId =threadIdx / WARP_SIZE;
            while( true )
            {
                int atomId;
                int cellId;
                if(atomIdId <level2PreciseWaterCnt)
                {
                    atomId = atomIdBuffer[atomIdBufferIdx + atomIdId];

                    cellId = atomId >> 16;
                    atomId = atomId & 0xffff;
                    level3scratch[threadIdx / WARP_SIZE] = atomId;

                    float4 a = atoms[atomId + localConstants.atomsOffset];
                    float sum = 0;

                    float dxa = localCell.x + level2root.w*(cellId & level2cellsperdimmask) - a.x;
                    float dya = localCell.y + level2root.w*((cellId>>(localConstants.maincelllevel + 2)) & level2cellsperdimmask) - a.y;
                    float dza = level2root.z +level2root.w*(cellId>>(4+ 2*localConstants.maincelllevel))- a.z;

                    float dy1 = dya*dya + dxa*dxa;
                    float dy2;
                    dy2 = (dya + 0.5f*level2root.w)*(dya + 0.5f*level2root.w) + dxa*dxa;
                    for(int zz = 0; zz < 8; zz++, dza += level2rootw125)
                    {
                        float d = dy1 + dza*dza;
                        sum += level2waterW3 / (d*d*d);
                        d = dy2 + dza*dza;
                        sum += level2waterW3 / (d*d*d);
                    }
                    sum /= 512.0f;

                    int scratchthread = threadIdx >> 2;
                    int scratchmask = threadIdx & 0x3;
                    if(scratchmask == 3)
                    {
                        scratch[scratchthread] = sum;
                    }
                    mem_fence(CLK_LOCAL_MEM_FENCE);
                    if(scratchmask == 1)
                    {
                        scratch[scratchthread] += sum;
                    }
                    mem_fence(CLK_LOCAL_MEM_FENCE);
                    if(scratchmask == 2)
                    {
                        scratch[scratchthread] += sum;
                    }
                    mem_fence(CLK_LOCAL_MEM_FENCE);
                    if(scratchmask == 0)
                    {
                        scratchmask = threadIdx & (WARP_SIZE - 1);
                        scratch[scratchthread] += sum;
                        mem_fence(CLK_LOCAL_MEM_FENCE);
                        if(scratchmask < 16)
                        {
                            scratch[scratchthread] += scratch[scratchthread + 4];
                        }
                        mem_fence(CLK_LOCAL_MEM_FENCE);
                        if(scratchmask < 8)
                        {
                            scratch[scratchthread] += scratch[scratchthread + 2];
                        }
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE);

                if(threadIdx == 0)
                {
                    for(int commitIdId = 0 ; commitIdId < GPU_BLOCK_SIZE/WARP_SIZE;commitIdId++) //switch to handling these additions in orderred bulk - save up to 10ms
                    {
                        if(commitIdId + atomIdId >= level2PreciseWaterCnt) break;

                        result[level3scratch[commitIdId] + localConstants.atomsOffset*integrationBlockCount + localConstants.atomPitch*integrationBlockId] += scratch[commitIdId*(WARP_SIZE/4)] + scratch[commitIdId*(WARP_SIZE/4) + 1];
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

                atomIdId += GPU_BLOCK_SIZE/WARP_SIZE;
                if(atomIdId /(GPU_BLOCK_SIZE/WARP_SIZE) > (level2PreciseWaterCnt - 1) /(GPU_BLOCK_SIZE/WARP_SIZE)) break;
            }
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE| CLK_GLOBAL_MEM_FENCE);
}


__kernel void
PowerBornKernel2(__global const float4* atoms, 
        __global const float4* water,
        __global const pGPUConstants* constants,
        __global float* result,
        __global  int* stack,
        __global int* atomIdBuffer,
        const int integrationBlocksScale
)
{
    __local int ss[16];
    __local int vv[16];
    __local int render[LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE];
    __local float4 level3[LEVEL3_SIZE];
    __local int level3scratch[LEVEL3_SIZE*4 + 1];

    __local float cellOffset  ;
    __local unsigned int level2cellsperdim ;
    __local unsigned int level2cellsperdim2 ;
    __local float4 level2root;

    pGPUConstants localConstants;

    __local float threshold;
    __local int level1CellDataOffset;
    __local int level1CellDataTotal;
    __local int level1CellId;
    __local int level1AtomIdId;
    __local int level1AtomIdCnt;
    __local int level1WaterIdId;
    __local int level1WaterIdCnt;
    __local int level1IntegIdCnt;
    float lastLevelW;
    float lastLevelW3;
    int threadIdx = get_local_id(0); // this is 20mus faster than calling get_local_id(0) all the time!
    // read init data
    int blockIdx = get_group_id(0);
    int integrationGroupId = get_group_id(1);
    int integrationBlockCount = get_num_groups(1);
    localConstants = constants[blockIdx];
    if(threadIdx == WARP_SIZEL)
    {
        unsigned int level1cellsperdim ;
        unsigned int level1cellsperdim2 ;
        level1cellsperdim = 1 << localConstants.maincelllevel;
        cellOffset  =  - ((float)level1cellsperdim / 2);
        threshold = 0.5f*localConstants.root.w / level1cellsperdim;
        threshold*= threshold*localConstants.integration_factor;
        level1cellsperdim2 = level1cellsperdim *level1cellsperdim;
        level2root.w = localConstants.root.w / level1cellsperdim;
        level2root.xyz = localConstants.root.xyz + (cellOffset)*level2root.www;
        cellOffset += 0.5f;

        level2cellsperdim = level1cellsperdim;
        level2cellsperdim2 = level1cellsperdim2;
        level1WaterIdCnt = 0;

        level1CellDataTotal = stack[localConstants.stackOffset - MAX_GPU_MAINCELLCOUNT*4];
        if(localConstants.maxLevel > 8)
        {
            level2cellsperdim <<=  2;
            level2cellsperdim2 <<= 4;
            level2root.w *= 0.25f;

            level1CellDataOffset = stack[localConstants.stackOffset - MAX_GPU_MAINCELLCOUNT*4 + 1];
            threshold /= 256.0f;
            cellOffset  = 0.5f - ((float)level2cellsperdim / 2);
        }
        else
        {
            level1CellDataOffset = MAX_GPU_MAINCELLCOUNT;
            threshold /= 16.0f;
        }

        int cellNum = 0;
        while ((cellNum < integrationGroupId) && (level1CellDataOffset < level1CellDataTotal))
        {
            level1CellDataOffset = stack[level1CellDataOffset + 2 + localConstants.stackOffset];
            cellNum++;
        }

    }
    // level1 generation lookups
    if(threadIdx < 16)
    {
        ss[threadIdx] = (threadIdx & 0x1) + ((threadIdx & 0x2) > 0) + ((threadIdx & 0x4) > 0) + ((threadIdx & 0x8) > 0);
        vv[threadIdx] = ((threadIdx>>1) & 0x3) + 3*((threadIdx>>3) & 0x1);
    }

    barrier(CLK_LOCAL_MEM_FENCE); 

    lastLevelW = level2root.w/32.0f;
    lastLevelW3 = lastLevelW*lastLevelW*lastLevelW;
    int atomIdBufferIdx = (get_group_id(0)*integrationBlockCount + integrationGroupId)*localConstants.atomPitch*integrationBlocksScale;
    int integrationBlockOffset = localConstants.atomPitch*integrationGroupId + localConstants.atomsOffset*integrationBlockCount;
    while (level1CellDataOffset < level1CellDataTotal)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if(threadIdx == 0)
        {
            level1AtomIdId = stack[level1CellDataOffset + localConstants.stackOffset]; // can do it as 1 16-byte read.
            level1WaterIdId = stack[level1CellDataOffset + 1 + localConstants.stackOffset];
            level1WaterIdCnt = stack[level1CellDataOffset + 2 + localConstants.stackOffset];
            level1CellId = stack[level1CellDataOffset + 3 + localConstants.stackOffset];
            level1AtomIdCnt = level1WaterIdId;
            level1CellDataOffset = level1WaterIdCnt; //safe to do this
            int cc = integrationBlockCount - 1;
            while ((cc > 0) && (level1CellDataOffset < level1CellDataTotal))
            {
                level1CellDataOffset = stack[level1CellDataOffset + 2 + localConstants.stackOffset];
                cc--;
            }

            int px = level1CellId%level2cellsperdim;
            int py = (level1CellId/level2cellsperdim)%level2cellsperdim;
            int pz = level1CellId/level2cellsperdim2;
            float lastlevelgrid = ((float)((LAST_LEVEL_GRID_SIZE-1)*0.5f ) )*lastLevelW;
            level2root.x = localConstants.root.x + (cellOffset + (float)px)*level2root.w - lastlevelgrid;
            level2root.y = localConstants.root.y + (cellOffset + (float)py)*level2root.w - lastlevelgrid;
            level2root.z = localConstants.root.z + (cellOffset + (float)pz)*level2root.w - lastlevelgrid;
            level1IntegIdCnt = 0;
        }
        // render 
        if(threadIdx >= WARP_SIZEL)
        {
            for(int atomid = threadIdx - WARP_SIZEL; atomid < LAST_LEVEL_GRID_SIZE*LAST_LEVEL_GRID_SIZE; atomid += GPU_BLOCK_SIZE_KERNEL2 - WARP_SIZEL) 
                render[atomid] = 0xffffffff;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        // renders atoms and water into render array
        // in/out render - each bit is set if corresponding last level cell contains water. starts with all water, renders atoms, then renders water
        // out level3 - xyzw coords of 64 cells that consist of 512 last level cells - this array used for coarse integrate
        // modifies level3scratch - for temp storage
        // modifies level1AtomIdId, level1WaterIdId - counts them up to level1AtomIdCnt, level1WaterIdCnt

        RenderCellAtomsAndWater(
                render, level3, level3scratch, &level1AtomIdId,
                level1WaterIdId, level1AtomIdCnt, level1WaterIdCnt, level2root, lastLevelW, lastLevelW3, threadIdx, localConstants, stack, atoms, water, ss, vv
        );


        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
        // computes coarse integrals for all atoms over level3 cells and marks cells that need precise integrals to be consumed by next function
        // out level1IntegIdCnt contains count of precise atom/cell pairs
        // out level3scratch contains first LEVEL3_SIZE*4 pairs
        // out stack contains the rest of pairs starting at stackPos, if any
        // in/out tempIntegrals adds coarse integral values to total integrals for atoms

        ComputeCoarseIntegral(
                &level1IntegIdCnt, level3scratch, atomIdBuffer, result,
                threadIdx, threshold, atomIdBufferIdx, integrationBlockOffset, localConstants, atoms, level3
        );

        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);


        // now do precise integrals for pairs selected by previous function
        // in/out tempIntegrals adds precise integral values for selected atom/cell pairs to total integrals for atoms
        // modifies scratch - for temp storage
        if(level1IntegIdCnt > 0)
        {
            /*
// this needs to be fixed to match code below
            ComputePreciseIntegral(
                result, scratch,
                level3scratch, atomIdBuffer, level1IntegIdCnt, threadIdx, atomIdBufferIdx, integrationBlockOffset, localConstants, level2root, lastLevelW, lastLevelW3, render, atoms
            );
             */
            // compiler fails at inlining... this is faster and uses less registers.
            {
                float3 localCell = level2root.xyz;
                localCell.x += (threadIdx & 0x7) * lastLevelW;
                localCell.y += ((threadIdx>>3) & 0x7) * lastLevelW;
                int zcoord = (threadIdx>>4) & 0x4;
                localCell.z += (zcoord) * lastLevelW;
                int atomIdId =0;
                while( atomIdId <level1IntegIdCnt)
                {
                    int atomId;
                    int mask;
                    int mask2;
                    if(atomIdId >= LEVEL3_SIZE)
                    {
                        atomId = atomIdBuffer[atomIdBufferIdx + atomIdId*4];
                        mask = atomIdBuffer[atomIdBufferIdx + atomIdId*4 + 1];
                        mask2 = atomIdBuffer[atomIdBufferIdx + atomIdId*4 + 2];
                    }
                    else
                    {
                        atomId = level3scratch[atomIdId*4];
                        mask = level3scratch[atomIdId*4 + 1];
                        mask2 = level3scratch[atomIdId*4 + 2];
                    }
                    float4 a = atoms[atomId + localConstants.atomsOffset];
                    a.xyz = localCell.xyz - a.xyz;
                    float sum2 = 0;
                    int cellIdOffset = 0;
                    if(!mask)
                    {
                        mask = mask2;
                        mask2 = 0;
                        cellIdOffset = 32;
                    }

                    while (mask )
                    {
                        uint cellId = 31 - clz(mask);
                        mask -= 1 << cellId;
                        cellId+=cellIdOffset;
                        if(mask == 0)
                        {
                            mask = mask2;
                            mask2 = 0;
                            cellIdOffset = 32;
                        }
                        float3 da;
                        da.x = cellId & 0x3;
                        da.y = (cellId>>2) & 0x3;
                        da.z = cellId>>4;
                        da.xyz = a.xyz + 8.0f*lastLevelW*da.xyz;
                        da.x=da.x*da.x + da.y*da.y;
                        int maskx = (threadIdx & 0x7) + ((cellId & 0x3) << 3);
                        int masky = ((threadIdx>>3) & 0x7) + ((cellId & 0xc) << 1);
                        int rendermask = render[masky*LAST_LEVEL_GRID_SIZE + maskx] >> (((cellId>>1) & 0x18) + (zcoord));
                        for (int z =0; z < 4; z++, da.z+= lastLevelW, rendermask >>= 1)
                        {
                            if(rendermask & 0x1)
                            {
                                float d = da.x + da.z*da.z;
                                sum2 += 1.0f / (d*d*d);
                            }
                        }
                    }
                    sum2*= lastLevelW3;
                    // reduce
                    int s1 = threadIdx & 0x1;
                    int s2 = threadIdx >> 1;
                    if(s1 == 0)
                        level3[s2 ].w = sum2;
                    mem_fence(CLK_LOCAL_MEM_FENCE);
                    if(s1 == 1)
                        level3[s2 ].w += sum2;
                    barrier(CLK_LOCAL_MEM_FENCE);
                    int offset = 32;
                    if(threadIdx < offset)
                    {
                        float res = result[atomId + integrationBlockOffset];
                        while (offset > 1)
                        {
                            if(threadIdx < offset)
                            {
                                level3[threadIdx].w += level3[threadIdx + offset].w;
                            }
                            mem_fence(CLK_LOCAL_MEM_FENCE);
                            offset >>= 1;
                        }
                        result[atomId + integrationBlockOffset] = res + (level3[0].w + level3[1].w);
                    }
                    barrier(CLK_LOCAL_MEM_FENCE);
                    atomIdId++;
                }
                barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

            }



        }
    }

    int divIdx = threadIdx / 8;
    int modIdx = threadIdx % 8;
    for(unsigned int atomid = 16*integrationGroupId; atomid < localConstants.atomCnt; atomid += 16*get_num_groups(1))
    {
        float3 tmp = atoms[atomid + divIdx + localConstants.atomsOffset].xyz - localConstants.root.xyz;
        float3 xp = 0.5f * localConstants.root.www - tmp;
        float3 xm = -0.5f * localConstants.root.www - tmp;
        float3 xpa;
        xpa.xyz = fabs(xp.xyz);
        float3 xma;
        xma.xyz = fabs(xm.xyz);
        float integral0 = 0;

        if(modIdx == 0)
        {
            float invr1=invR(xma.x,xpa.y,xpa.z);
            float angley0 = asin(invr1* xpa.y);
            float anglez0 = asin(invr1* xpa.z);
            integral0 += squareIntegral(xm.x, xp.y, xp.z, angley0, anglez0); // back
            float angley14 = asin(invr1* xma.x);
            integral0 += squareIntegral(xp.y, xm.x, xp.z, angley14, anglez0);
            integral0 += squareIntegral(xp.z, xp.y, xm.x, angley0, angley14);
            level3[divIdx].x = integral0;
        }
        else if(modIdx == 1)
        {
            float invr2=invR(xma.x,xpa.y,xma.z);
            float angley1 = asin(invr2* xpa.y);
            float anglez1 = asin(invr2* xma.z);
            integral0 -= squareIntegral(xm.x, xp.y, xm.z, angley1, anglez1);
            float angley15 = asin(invr2* xma.x);
            integral0 -= squareIntegral(xp.y, xm.x, xm.z, angley15, anglez1);
            integral0 -= squareIntegral(xm.z, xp.y, xm.x, angley1, angley15);
            level3[divIdx].y = integral0;
        }
        else if(modIdx == 2)
        {
            float invr3=invR(xma.x,xma.y,xpa.z);
            float angley2 = asin(invr3* xma.y);
            float anglez2 = asin(invr3* xpa.z);
            integral0 -= squareIntegral(xm.x, xm.y, xp.z, angley2, anglez2);
            float angley10 = asin(invr3* xma.x);
            integral0 -= squareIntegral(xm.y, xm.x, xp.z, angley10, anglez2);
            integral0 -= squareIntegral(xp.z, xm.y, xm.x, angley2, angley10);
            level3[divIdx].z = integral0;
        }
        else if(modIdx == 3)
        {
            float invr4=invR(xma.x,xma.y,xma.z);
            float angley3 = asin(invr4* xma.y);
            float anglez3 = asin(invr4* xma.z);
            integral0 += squareIntegral(xm.x, xm.y, xm.z, angley3, anglez3);
            float angley11 = asin(invr4* xma.x);
            integral0 += squareIntegral(xm.y, xm.x, xm.z, angley11, anglez3);
            integral0 += squareIntegral(xm.z, xm.y, xm.x, angley3, angley11);
            level3[divIdx].w = integral0;
        }
        else if(modIdx == 4)
        {
            float invr5=invR(xpa.x,xpa.y,xpa.z);
            float angley4 = asin(invr5* xpa.y);
            float anglez4 = asin(invr5* xpa.z);
            integral0 -= squareIntegral(xp.x, xp.y, xp.z, angley4, anglez4); //front
            float angley12 = asin(invr5* xpa.x);
            integral0 -= squareIntegral(xp.y, xp.x, xp.z, angley12, anglez4); //right
            integral0 -= squareIntegral(xp.z, xp.y, xp.x, angley4, angley12); //top
            level3[32 + divIdx].x = integral0;
        }
        else if(modIdx == 5)
        {
            float invr6=invR(xpa.x,xpa.y,xma.z);
            float angley5 = asin(invr6* xpa.y);
            float anglez5 = asin(invr6* xma.z);
            integral0 += squareIntegral(xp.x, xp.y, xm.z, angley5, anglez5);
            float angley13 = asin(invr6* xpa.x);
            integral0 += squareIntegral(xp.y, xp.x, xm.z, angley13, anglez5);
            integral0 += squareIntegral(xm.z, xp.y, xp.x, angley5, angley13); //bottom
            level3[32 + divIdx].y = integral0;
        }
        else if(modIdx == 6)
        {
            float invr7=invR(xpa.x,xma.y,xpa.z);
            float angley6 = asin(invr7* xma.y);
            float anglez6 = asin(invr7* xpa.z);
            integral0 += squareIntegral(xp.x, xm.y, xp.z, angley6, anglez6);
            float angley8 = asin(invr7* xpa.x);
            integral0 += squareIntegral(xm.y, xp.x, xp.z, angley8, anglez6); //left
            integral0 += squareIntegral(xp.z, xm.y, xp.x, angley6, angley8);
            level3[32 + divIdx].z = integral0;
        }
        else if(modIdx == 7)
        {
            float invr8=invR(xpa.x,xma.y,xma.z);
            float angley7 = asin(invr8* xma.y);
            float anglez7 = asin(invr8* xma.z);
            integral0 -= squareIntegral(xp.x, xm.y, xm.z, angley7, anglez7);
            float angley9 = asin(invr8* xpa.x);
            integral0 -= squareIntegral(xm.y, xp.x, xm.z, angley9, anglez7);
            integral0 -= squareIntegral(xm.z, xm.y, xp.x, angley7, angley9);
            level3[32 + divIdx].w = integral0;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        if(threadIdx < 16)
        {
            float br = result[atomid + threadIdx + integrationBlockOffset]
                              + ( level3[threadIdx].x
                                      + level3[threadIdx].y
                                      + level3[threadIdx].z
                                      + level3[threadIdx].w
                                      + level3[32 + threadIdx].x
                                      + level3[32 + threadIdx].y
                                      + level3[32 + threadIdx].z
                                      + level3[32 + threadIdx].w);

            result[atomid + threadIdx + integrationBlockOffset] = br;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    barrier(CLK_LOCAL_MEM_FENCE| CLK_GLOBAL_MEM_FENCE);
}


__kernel void
PowerBornKernelAggregate(  
        __global const float4* atoms,
        __global const pGPUConstants* constants,
        __global float* result,
        const int integrationBlockCount
)
{
    int threadIdx = get_local_id(0); // this is 20mus faster than calling get_local_id(0) all the time!
    int blockIdx = get_group_id(0);
    int atomCnt = constants[blockIdx].atomCnt;
    int atomOffset = constants[blockIdx].atomsOffset;
    int atomId = threadIdx + get_group_id(1)*GPU_BLOCK_SIZE;
    int atomPitch = constants[blockIdx].atomPitch;
    if(atomId < atomCnt)
    {
        float br = 0;
        for(int i = 0; i < integrationBlockCount; ++i) br += result[atomOffset*integrationBlockCount + atomId + i*atomPitch];

        br = pow(3.0f / 4.0f / M_PI * br, 0.3333333333333333f);
        br = (1.0f) / ( FIT_PARAM_A * br + FIT_PARAM_B);
        float rad = atoms[atomId + atomOffset].w - 1.4f;
        br = max(br, rad);
        result[atomOffset*integrationBlockCount + atomId] = br;
    }
    barrier(CLK_LOCAL_MEM_FENCE| CLK_GLOBAL_MEM_FENCE);
}


__kernel void
WaterKernel(__global const float4* atoms, 
        __global pGPUConstants* constants,
        __global const float4* waterSeeds,
        __global int* waterMask,
        __global int2* buffersizes)
{
    __local float4 srcatoms[GPU_BLOCK_SIZE/WARP_SIZE];
    __local int activeCnt;
    __local float4 dstatoms[GPU_BLOCK_SIZE];
    __local int active[GPU_BLOCK_SIZE];
    __local unsigned int localwatermask[(WATERMASK_SIZE/32)*(GPU_BLOCK_SIZE/WARP_SIZE)];
    __local unsigned int waterCnt;
    __local unsigned int stackCnt;
    __local int waterSeedsCnt;
    __local float preciseres;
    __local float preciserate;
    __local float4 root;
    __local int atomLimit;
    __local int atomCnt;
    __local int atomsOffset;
    __local int atomid;

    int threadIdx = get_local_id(0);
    // initial population of coarse grid
    if(threadIdx == 0)
    {
        waterCnt = 0;
        stackCnt = 0;
        root = constants[get_group_id(0)].root;
        atomCnt = constants[get_group_id(0)].atomCnt;
        atomsOffset = constants[get_group_id(0)].atomsOffset;
        int levels = constants[get_group_id(0)].maxLevel;
        preciseres = 1.0f / (root.w / (1 << (levels - 5)));
        preciserate = (2.0f * 1.4f) * preciseres;
        atomLimit = (get_group_id(1)+1)*(GPU_BLOCK_SIZE/WARP_SIZE);
        if (atomLimit > atomCnt) atomLimit = atomCnt;
        atomid = get_group_id(1)*(GPU_BLOCK_SIZE/WARP_SIZE);
    }

    // reduce min
    barrier(CLK_LOCAL_MEM_FENCE);

    int localid = threadIdx%WARP_SIZE;
    int midid = threadIdx/WARP_SIZE; 
    int localWaterCnt = 0; 
    int localStackCnt = 0; 
    unsigned int shift = 1 << localid;
    if(threadIdx < GPU_BLOCK_SIZE/WARP_SIZE && (atomid +threadIdx < atomLimit))
    {
        srcatoms[threadIdx] = atoms[atomid + threadIdx + atomsOffset];
    }
    if(localid < WATERMASK_SIZE/32) localwatermask[localid + midid*(WATERMASK_SIZE/32)] = 0xffffffff;
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int testatomid = 0; testatomid < atomCnt; testatomid += GPU_BLOCK_SIZE)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        if(testatomid + threadIdx < atomCnt)
        {
            dstatoms[threadIdx] = atoms[testatomid + threadIdx + atomsOffset];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int activeatomid = 0;activeatomid < GPU_BLOCK_SIZE/WARP_SIZE;activeatomid ++)
        {
            if(atomid + activeatomid >= atomLimit) break;
            if(threadIdx == 0)
            {
                activeCnt = 0;
                if(srcatoms[activeatomid].w < WATERMASK_THRESHOLD)
                {
                    waterSeedsCnt = WATERMASK_SMALL;
                }
                else
                {
                    waterSeedsCnt = WATERMASK_LARGE;
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if(testatomid+threadIdx < atomCnt && testatomid+threadIdx != atomid + activeatomid)
            {
                float3 d = dstatoms[threadIdx].xyz - srcatoms[activeatomid].xyz;
                float dw = dstatoms[threadIdx].w + srcatoms[activeatomid].w;
                if(d.x*d.x + d.y*d.y+d.z*d.z < dw*dw)
                {
                    int pos = atomic_add(&activeCnt,1);
                    active[pos] = threadIdx;
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            int maskadr = (WATERMASK_SIZE/32)*activeatomid + (threadIdx/32);
            int wid = threadIdx;
            for(;wid < waterSeedsCnt;wid += GPU_BLOCK_SIZE,maskadr+=GPU_BLOCK_SIZE/32)
            {
                if(localwatermask[maskadr] & shift)
                {
                    float3 srcw;
                    if(srcatoms[activeatomid].w < WATERMASK_THRESHOLD)
                    {
                        srcw = waterSeeds[wid].xyz;
                    }
                    else
                    {
                        srcw = waterSeeds[wid + WATERMASK_SMALL].xyz;
                    }
                    srcw = srcatoms[activeatomid].xyz + srcatoms[activeatomid].www*srcw.xyz;
                    int testbufferid = 0;
                    int is_covered = 0;
                    while ((testbufferid < activeCnt) & !is_covered)
                    {
                        int dstid = active[testbufferid];
                        float3 dist_w = srcw - dstatoms[dstid].xyz;
                        is_covered |= dist_w.x*dist_w.x + dist_w.y*dist_w.y + dist_w.z*dist_w.z < dstatoms[dstid].w*dstatoms[dstid].w;
                        testbufferid++;
                    }
                    if(is_covered)
                    {
                        atomic_and(&localwatermask[maskadr], 0xffffffff - shift);
                    }
                    if(!is_covered && testatomid + GPU_BLOCK_SIZE >= atomCnt && testbufferid == activeCnt)
                    {
                        // this water survived filtering, compute stack allocation
                        srcw.xyz -= root.xyz;
                        srcw.x = (srcw.x + 1.4f) * preciseres;
                        srcw.y = (srcw.y + 1.4f) * preciseres;
                        srcw.z = (srcw.z + 1.4f) * preciseres;
                        int c = (srcw.x - floor(srcw.x) <= preciserate) + (srcw.y - floor(srcw.y)<= preciserate) + (srcw.z - floor(srcw.z) <= preciserate);
                        localStackCnt += 1 << c;
                        localWaterCnt++;
                    }
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }
    mem_fence(CLK_LOCAL_MEM_FENCE);
    if(localid < WATERMASK_SIZE/32 && (atomid + midid < atomLimit))
    {
        int a = atomid + midid + atomsOffset;
        waterMask[a*(WATERMASK_SIZE/32) + localid] = localwatermask[localid + midid*(WATERMASK_SIZE/32)];
    }

    atomic_add(&waterCnt,localWaterCnt);
    atomic_add(&stackCnt,localStackCnt);
    barrier(CLK_LOCAL_MEM_FENCE);
    if(threadIdx == 0)
    {
        int pos = get_group_id(0)*get_global_size(1) + get_group_id(1);
        buffersizes[pos].x = waterCnt;
        buffersizes[pos].y = stackCnt;
    }
    barrier(CLK_GLOBAL_MEM_FENCE);
}

__kernel void
EnergyKernel(__global const float4* atoms,
        __global const pGPUEnergyConstants* energyconstants,
        __global const pGPUConstants* constants,
        const int integrationBlockCount,
        __global const float* result,
        __global const float* resultArea,
        __global float4* resultEnergy,
        __global const float* eps,
        __global const float* sig,
        __global const float* charge,
        __global const unsigned int* excludedmask,
        __global const int4* dihdefs,
        __global const float4* dihparams_v,
        __global const float4* dihparams_t,
        __global const float4* dihparams_m,
        __global int2* buffersizes)
{
    __local float4 params[GPU_BLOCK_SIZE];
    __local float3 params2[16];

    const unsigned int threadIdx = get_local_id(0);
    const int blockIdx = get_group_id(0);
    const unsigned int atomOffset = constants[blockIdx].atomsOffset;
    const unsigned int atomCnt = constants[blockIdx].atomCnt;
    if(get_group_id(1) > 0)
    { 
        const unsigned int atomPitch = constants[blockIdx].atomPitch;
        // begin LJ, Coulomb, GB
        float excluded_param_co = energyconstants->excluded_param_co;
        float excluded_param_lj = energyconstants->excluded_param_lj;

        float coulomb = 0;
        float lj = 0;
        float gb = 0;
        unsigned int i = (get_group_id(1) - 1)*GPU_BLOCK_SIZE;
        {
            float3 a = {0.0f, 0.0f, 0.0f};
            float br1, charge1, eps1, sig1;
            br1 = sig1 = 1.0f;
            charge1 = eps1 = 0.0f;

            if(i + threadIdx < atomCnt)
            {
                const unsigned int atomid1 = i + threadIdx ;
                a = atoms[atomid1 + atomOffset].xyz;
                br1 = result[atomid1 + atomOffset*integrationBlockCount];
                charge1 = charge[atomid1];
                eps1 = eps[atomid1];
                sig1 = sig[atomid1];
            }
            for(unsigned int block = 0; block<atomCnt; block+=16)
            {
                if(threadIdx < 16)
                {
                    if(block + threadIdx < atomCnt)
                    {
                        params[threadIdx].xyz = atoms[threadIdx + block + atomOffset].xyz;
                        params[threadIdx].w = result[threadIdx + block + atomOffset*integrationBlockCount];
                        params2[threadIdx].x = charge[threadIdx + block];
                        params2[threadIdx].y = eps[threadIdx + block];
                        params2[threadIdx].z = sig[threadIdx + block];
                    }
                    else 
                    {
                        params[threadIdx].xyzw = 1.0f;
                        params2[threadIdx].xyz = 0.0f;
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE);

                if(i + threadIdx < atomCnt)
                {
                    unsigned int mymask = excludedmask[((i + threadIdx) * atomPitch + block) / 16];
                    for(unsigned int j=0; j<16; j++)
                    {
                        unsigned int mask = mymask;
                        mask >>= (j * 2);

                        float charge12 = params2[j].x * charge1;
                        float eps12 = params2[j].y * eps1;
                        float sig12 = params2[j].z + sig1;

                        float3 b = params[j].xyz - a.xyz;
                        float br12 = params[j].w * br1;

                        float rsqr = ((b.x * b.x + b.y * b.y) + (b.z * b.z + 0.000001f));
                        float inv_rsqr = 1.0f / rsqr;
                        charge12 *= sqrt(inv_rsqr);
                        float r6 = sig12 * sig12 * inv_rsqr;
                        r6 *= r6 * r6;
                        float lj_res = eps12 * (r6 * r6 - r6);

                        lj_res = mask & 1? 0.0f : lj_res;
                        const float ex_param = mask & 1 ? 0.0f : 1.0f;
                        const float ex_coulomb_param = mask & 2 ? excluded_param_co : ex_param;
                        const float ex_lj_param = mask & 2 ? excluded_param_lj : 1.0f;

                        coulomb += charge12 * ex_coulomb_param;
                        lj += lj_res * ex_lj_param;
                        gb += charge12 / sqrt(1.0f + br12 * inv_rsqr * exp(-rsqr/(4.0f * br12)));
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }
        }
        params[threadIdx].x = coulomb;
        params[threadIdx].y = lj;
        params[threadIdx].z = gb;
        params[threadIdx].w = 0.0f;
        if(i + threadIdx < atomCnt)
        {
            params[threadIdx].w = resultArea[atomOffset + i + threadIdx];
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // reduce
        int offset = 256 / 2;
        while(offset)
        {
            if(threadIdx <offset)
            {
                params[threadIdx] += params[threadIdx + offset];
            }
            offset /= 2;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        if(threadIdx == 0)
        {
            float4 p = params[0];
            p.x *= energyconstants->coulomb_param;
            p.z *= energyconstants->gb_param;
            p.w *= energyconstants->npse_param;
            p.y *= (4.0f * 0.5f);
            unsigned int store = blockIdx;
            store *= 1 + (atomCnt + GPU_BLOCK_SIZE - 1)/GPU_BLOCK_SIZE;
            store += get_group_id(1);
            resultEnergy[store] = p;
            int pos = get_group_id(0) * get_global_size(1) + get_group_id(1);
            buffersizes[pos] = 1;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    else
    {
        // Dihedral terms
        float4 dih_result = {0.0f, 0.0f, 0.0f, 0.0f};
        unsigned int dihCnt = energyconstants->dihCnt;
        for(int dihId = threadIdx; dihId < dihCnt; dihId += GPU_BLOCK_SIZE)  
        {
            int4 dihIds = dihdefs[dihId];
            float3 c1 = atoms[dihIds.x + atomOffset].xyz;
            float3 c2 = atoms[dihIds.y + atomOffset].xyz;
            float3 c3 = atoms[dihIds.z + atomOffset].xyz;
            float3 c4 = atoms[dihIds.w + atomOffset].xyz;
            float3 b1,b2,b3;
            b1 = c2-c1;
            b2 = c3-c2;
            b3 = c4-c3;
            float4 theta;
            theta.xyzw = atan2(length(b2)*dot(b1,cross(b2,b3)), dot(cross(b1,b2),cross(b2, b3)));
            float4 cosine = dihparams_m[dihId] * theta - dihparams_t[dihId];
            cosine.x = cos(cosine.x);
            cosine.y = cos(cosine.y);
            cosine.z = cos(cosine.z);
            cosine.w = cos(cosine.w);
            dih_result += dihparams_v[dihId] * (1.0f + cosine);
        }
        params[threadIdx].x = dih_result.x + dih_result.y + dih_result.z + dih_result.w;
        barrier(CLK_LOCAL_MEM_FENCE);

        // reduce
        int offset = 256 / 2;
        while(offset)
        {
            if(threadIdx <offset)
            {
                params[threadIdx].x += params[threadIdx + offset].x;
            }
            offset /= 2;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        if(threadIdx == 0)
        {
            float4 p = params[0];
            p.yzw = 0;
            unsigned int store = blockIdx;
            store *= 1 + (atomCnt + GPU_BLOCK_SIZE - 1)/GPU_BLOCK_SIZE;
            resultEnergy[store] = p;
            int pos = get_group_id(0) * get_global_size(1) + get_group_id(1);
            buffersizes[pos].x = 1;
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE| CLK_GLOBAL_MEM_FENCE);
}
