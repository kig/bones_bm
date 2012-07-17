/*

Benchmark that applies bones to a mesh, uses SSE vectorization and
OpenMP for parallelism.

g++ -o bones -O3 -lrt -msse3 -fopenmp bones.cpp # compile OpenMP version
g++ -o bones_ST -O3 -lrt -msse3 bones.cpp # compile single-threaded version
./bones # run SSE version
./bones scalar # run scalar version

*/

#include "sse.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define VERT_COUNT 125000
#define BONES_COUNT 200

using namespace std;

typedef unsigned long long uint64_t;

struct weight_t {
  float weights[2];
  int bones[2];
  int length;
} typedef weight;

weight randomWeights(int boneCount) {
  weight w;
  w.length = 1 + (rand() / (RAND_MAX/2));
  for (int i=0; i<w.length; i++) {
    w.bones[i] = (rand() / (RAND_MAX/boneCount));
    w.weights[i] = ((float)rand()) / RAND_MAX;
  }
  return w;
}

float4 mul_m4x4_v4(float4 *matrix, float4 v) {
  float4 dst = matrix[0] * float4(v.x());
  dst += matrix[1] * float4(v.y());
  dst += matrix[2] * float4(v.z());
  dst += matrix[3] * float4(v.w());
  return dst;
}

/**
  @param {Array<Float32Array>} dstVertices Vertex array to contain the transformed vertices
  @param {Array<Float32Array>} srcVertices Vertex array with untransformed vertices
  @param {Array<Object<bones: Uint32Array, weights: Float32Array>>} weights
    Array that maps vertices to a list of (bone, weight) pairs.
  @param {Array<Array<Float32Array>>} bones Bone matrices.
  */
void applyBones(float4 *dstVertices, float4 *srcVertices, weight* weights, float4 **bones) {
  #pragma omp parallel for
  for (int i = 0; i < VERT_COUNT; i++) {
    float4 v = srcVertices[i];
    if (weights[i].length == 0) {
      dstVertices[i] = v;
    } else {
      /* Sum up the transformed verts blended with their weights,
         normalize to 1.0 total weight */
      float4 dv(0.0f);
      weight w = weights[i];
      float totalWeight = 0.0f;
      for (int j = 0; j < w.length; j++) {
        float4 *bone = bones[w.bones[j]];
        float4 tmp = mul_m4x4_v4(bone, v);
        // Scale transformed vert by its weight
        tmp *= float4(w.weights[j]);
        dv += tmp;
        totalWeight += w.weights[j];
      }
      dv /= float4(totalWeight);
      dstVertices[i] = dv;
    }
  }
};

void applyBones_scalar(float *dstVertices, float *srcVertices, weight* weights, float **bones) {
  #pragma omp parallel for
  for (int i = 0; i < VERT_COUNT; i++) {
    int off = i * 4;
    if (weights[i].length == 0) {
      dstVertices[off+0] = srcVertices[off+0];
      dstVertices[off+1] = srcVertices[off+1];
      dstVertices[off+2] = srcVertices[off+2];
      dstVertices[off+3] = srcVertices[off+3];
    } else {
      /* Sum up the transformed verts blended with their weights,
         normalize to 1.0 total weight */
      dstVertices[off+0] = 0.0f;
      dstVertices[off+1] = 0.0f;
      dstVertices[off+2] = 0.0f;
      dstVertices[off+3] = 0.0f;
      weight w = weights[i];
      float totalWeight = 0.0f;
      for (int j = 0; j < w.length; j++) {
        float *bone = bones[w.bones[j]];
        float wt = w.weights[j];
        totalWeight += wt;
        float x = srcVertices[off+0], y = srcVertices[off+1], z = srcVertices[off+2], w = srcVertices[off+3];
        dstVertices[off+0] += wt * (bone[0] * x + bone[4] * y + bone[8] * z + bone[12] * w);
        dstVertices[off+1] += wt * (bone[1] * x + bone[5] * y + bone[9] * z + bone[13] * w);
        dstVertices[off+2] += wt * (bone[2] * x + bone[6] * y + bone[10] * z + bone[14] * w);
        dstVertices[off+3] += wt * (bone[3] * x + bone[7] * y + bone[11] * z + bone[15] * w);
      }
      dstVertices[off+0] /= totalWeight;
      dstVertices[off+1] /= totalWeight;
      dstVertices[off+2] /= totalWeight;
      dstVertices[off+3] /= totalWeight;
    }
  }
};

int main (int argc, char* argv[]) {
  timespec t0, t1;
  float4 *srcVerts = (float4*)malloc(VERT_COUNT*sizeof(float4));
  float4 *dstVerts = (float4*)malloc(VERT_COUNT*sizeof(float4));
  weight *weights = (weight*)malloc(VERT_COUNT*sizeof(weight));

  for (int i = 0; i < VERT_COUNT; i++) {
    srcVerts[i] = float4(1.0);
    dstVerts[i] = float4(2.0);
    weights[i] = randomWeights(BONES_COUNT);
  }

  float4 *bonesBuffer = (float4*)malloc(BONES_COUNT*4*sizeof(float4));
  float4 **bones = (float4**)malloc(BONES_COUNT*sizeof(float4*));
  for (int i = 0; i < BONES_COUNT; i++) {
    bones[i] = bonesBuffer + (i*4);
    bones[i][0] = float4(3.0);
    bones[i][1] = float4(3.0);
    bones[i][2] = float4(3.0);
    bones[i][3] = float4(3.0);
  }

  if (argc > 1) {
    printf("Running scalar version...\n");
  } else {
    printf("Running SSE version...\n");
  }

  clock_gettime(CLOCK_REALTIME, &t0);
  int i;
  /* Main benchmark loop */
  for (i = 0; i < 100; i++) {
    if (argc > 1) {
      applyBones_scalar((float*)dstVerts, (float*)srcVerts, weights, (float**)bones);
    } else {
      applyBones(dstVerts, srcVerts, weights, bones);
    }
    for (int j = 0; j < BONES_COUNT; j++) {
      bones[j][0] = float4((float)i);
    }
  }
  clock_gettime(CLOCK_REALTIME, &t1);

  uint64_t elapsed = ((uint64_t)t1.tv_sec*1000000LL + (uint64_t)t1.tv_nsec/1000LL) -
                     ((uint64_t)t0.tv_sec*1000000LL + (uint64_t)t0.tv_nsec/1000LL);

  printf("elapsed: %llu us\n", elapsed / i);

  float4 v = dstVerts[0];
  printf("dstVerts[0]: %f, %f, %f, %f\n", v.x(), v.y(), v.z(), v.w());

  free(weights);
  free(srcVerts);
  free(dstVerts);
  free(bones);
  free(bonesBuffer);

  return 0;
}