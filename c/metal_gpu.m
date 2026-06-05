#include "gtc.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <stdlib.h>
#include <string.h>

typedef struct {
  int mi;
  int mpsi;
  int mzeta;
  int single_zeta_cell;
  float delr;
  float delz;
  float smu_inv;
  float pi2_inv;
  float zetamin;
  float a0;
  float a1;
} GtcMetalChargeParams;

typedef struct {
  int mi;
  int mimax;
  int mpsi;
  int mzeta;
  int irk;
  int linear_orbit;
  int has_vdrtmp;
  int write_wpi;
  float delr;
  float pi2;
  float psimax;
  float cmratio;
  float cinv;
  float vthi;
  float ainv;
  float sbound;
  float dtime;
  float a;
  float a0;
  float q0;
  float q1;
  float q2;
  float qion;
  float aion;
  float kappati;
  float kappan;
  float paranl;
  float nonlinear;
  float flow0;
  float rc;
  float rw;
  float gyroradius;
} GtcMetalPushiParams;

static const char *gtc_metal_source =
"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"enum {\n"
"  GTC_Z_PSI = 0,\n"
"  GTC_Z_THETA = 1,\n"
"  GTC_Z_ZETA = 2,\n"
"  GTC_Z_U = 3,\n"
"  GTC_Z_WEIGHT = 4,\n"
"  GTC_Z_MU = 5,\n"
"  GTC_NPARAM = 6\n"
"};\n"
"\n"
"struct ChargeParams {\n"
"  int mi;\n"
"  int mpsi;\n"
"  int mzeta;\n"
"  int single_zeta_cell;\n"
"  float delr;\n"
"  float delz;\n"
"  float smu_inv;\n"
"  float pi2_inv;\n"
"  float zetamin;\n"
"  float a0;\n"
"  float a1;\n"
"};\n"
"\n"
"struct PushiParams {\n"
"  int mi;\n"
"  int mimax;\n"
"  int mpsi;\n"
"  int mzeta;\n"
"  int irk;\n"
"  int linear_orbit;\n"
"  int has_vdrtmp;\n"
"  int write_wpi;\n"
"  float delr;\n"
"  float pi2;\n"
"  float psimax;\n"
"  float cmratio;\n"
"  float cinv;\n"
"  float vthi;\n"
"  float ainv;\n"
"  float sbound;\n"
"  float dtime;\n"
"  float a;\n"
"  float a0;\n"
"  float q0;\n"
"  float q1;\n"
"  float q2;\n"
"  float qion;\n"
"  float aion;\n"
"  float kappati;\n"
"  float kappan;\n"
"  float paranl;\n"
"  float nonlinear;\n"
"  float flow0;\n"
"  float rc;\n"
"  float rw;\n"
"  float gyroradius;\n"
"};\n"
"\n"
"static inline int clamp_int(int value, int lo, int hi) {\n"
"  return min(max(value, lo), hi);\n"
"}\n"
"\n"
"static inline float modulo_positive(float x, float period) {\n"
"  if (x >= period) {\n"
"    x -= period;\n"
"    if (x >= period) x = fmod(x, period);\n"
"  } else if (x < 0.0f) {\n"
"    x += period;\n"
"    if (x < 0.0f) {\n"
"      x = fmod(x, period);\n"
"      if (x < 0.0f) x += period;\n"
"    }\n"
"  }\n"
"  return x;\n"
"}\n"
"\n"
"static inline uint particle_index(uint param, uint m) {\n"
"  return m * (uint)GTC_NPARAM + param;\n"
"}\n"
"\n"
"static inline uint gyro_index(uint larmor, uint ij) {\n"
"  return ij * 4u + larmor;\n"
"}\n"
"\n"
"static inline uint larmor_particle_index(uint larmor, uint m) {\n"
"  return m * 4u + larmor;\n"
"}\n"
"\n"
"static inline uint evector_index(constant PushiParams &p, uint component, uint kz, uint ij) {\n"
"  return ((ij * (uint)(p.mzeta + 1) + kz) * 3u) + component;\n"
"}\n"
"\n"
"kernel void chargei_prepare_kernel(device const float *zion [[buffer(0)]],\n"
"                                  device const float *pgyro [[buffer(1)]],\n"
"                                  device const float *tgyro [[buffer(2)]],\n"
"                                  device const int *mtheta [[buffer(3)]],\n"
"                                  device const int *igrid [[buffer(4)]],\n"
"                                  device const float *qtinv [[buffer(5)]],\n"
"                                  device const float *delt [[buffer(6)]],\n"
"                                  device int *kzion [[buffer(7)]],\n"
"                                  device float *wzion [[buffer(8)]],\n"
"                                  device float *wpion [[buffer(9)]],\n"
"                                  device int *jtion0 [[buffer(10)]],\n"
"                                  device int *jtion1 [[buffer(11)]],\n"
"                                  device float *wtion0 [[buffer(12)]],\n"
"                                  device float *wtion1 [[buffer(13)]],\n"
"                                  constant ChargeParams &p [[buffer(14)]],\n"
"                                  uint m [[thread_position_in_grid]]) {\n"
"  if ((int)m >= p.mi) return;\n"
"  const device float *zp = zion + m * (uint)GTC_NPARAM;\n"
"  const float psitmp = zp[GTC_Z_PSI];\n"
"  const float thetatmp = zp[GTC_Z_THETA];\n"
"  const float zetatmp = zp[GTC_Z_ZETA];\n"
"  const float rhoi = zp[GTC_Z_MU] * p.smu_inv;\n"
"  const float r = sqrt(2.0f * psitmp);\n"
"  const int ip = clamp_int((int)((r - p.a0) * p.delr + 0.5f), 0, p.mpsi);\n"
"  const int jt = clamp_int((int)(thetatmp * p.pi2_inv * delt[ip] + 0.5f), 0, mtheta[ip]);\n"
"  const int ipjt = igrid[ip] + jt;\n"
"  const float wz1 = (zetatmp - p.zetamin) * p.delz;\n"
"  if (p.single_zeta_cell != 0) {\n"
"    kzion[m] = 0;\n"
"    wzion[m] = wz1;\n"
"  } else {\n"
"    const int kk = clamp_int((int)wz1, 0, p.mzeta - 1);\n"
"    kzion[m] = kk;\n"
"    wzion[m] = wz1 - (float)kk;\n"
"  }\n"
"  for (uint larmor = 0; larmor < 4u; larmor++) {\n"
"    const float rdum = p.delr * max(0.0f, min(p.a1 - p.a0,\n"
"        r + rhoi * pgyro[gyro_index(larmor, (uint)ipjt)] - p.a0));\n"
"    const int ii = clamp_int((int)rdum, 0, p.mpsi - 1);\n"
"    const float wp1 = rdum - (float)ii;\n"
"    const uint o = larmor_particle_index(larmor, m);\n"
"    wpion[o] = wp1;\n"
"\n"
"    const float tflr = thetatmp + rhoi * tgyro[gyro_index(larmor, (uint)ipjt)];\n"
"    float tdum = p.pi2_inv * (tflr - zetatmp * qtinv[ii]) + 10.0f;\n"
"    tdum = (tdum - (float)((int)tdum)) * delt[ii];\n"
"    const int j00 = clamp_int((int)tdum, 0, mtheta[ii] - 1);\n"
"    jtion0[o] = igrid[ii] + j00;\n"
"    wtion0[o] = tdum - (float)j00;\n"
"\n"
"    const int im = ii + 1;\n"
"    tdum = p.pi2_inv * (tflr - zetatmp * qtinv[im]) + 10.0f;\n"
"    tdum = (tdum - (float)((int)tdum)) * delt[im];\n"
"    const int j01 = clamp_int((int)tdum, 0, mtheta[im] - 1);\n"
"    jtion1[o] = igrid[im] + j01;\n"
"    wtion1[o] = tdum - (float)j01;\n"
"  }\n"
"}\n"
"\n"
"kernel void pushi_general_kernel(device float *zion [[buffer(0)]],\n"
"                                device float *zion0 [[buffer(1)]],\n"
"                                device float *wpi [[buffer(2)]],\n"
"                                device const float *evector [[buffer(3)]],\n"
"                                device const int *kzion [[buffer(4)]],\n"
"                                device const int *jtion0 [[buffer(5)]],\n"
"                                device const int *jtion1 [[buffer(6)]],\n"
"                                device const float *wzion [[buffer(7)]],\n"
"                                device const float *wpion [[buffer(8)]],\n"
"                                device const float *wtion0 [[buffer(9)]],\n"
"                                device const float *wtion1 [[buffer(10)]],\n"
"                                device const float *temp_inv [[buffer(11)]],\n"
"                                device const float *vdrtmp [[buffer(12)]],\n"
"                                constant PushiParams &p [[buffer(13)]],\n"
"                                uint m [[thread_position_in_grid]]) {\n"
"  if ((int)m >= p.mi) return;\n"
"\n"
"  float e1 = 0.0f;\n"
"  float e2 = 0.0f;\n"
"  float e3 = 0.0f;\n"
"  const uint kk = (uint)kzion[m];\n"
"  const float wz1 = wzion[m];\n"
"  const float wz0 = 1.0f - wz1;\n"
"  for (uint larmor = 0; larmor < 4u; larmor++) {\n"
"    const uint o = larmor_particle_index(larmor, m);\n"
"    uint ij = (uint)jtion0[o];\n"
"    const float wp0 = 1.0f - wpion[o];\n"
"    const float wt00 = 1.0f - wtion0[o];\n"
"    const float wt10 = 1.0f - wt00;\n"
"    float coeff = wp0 * wt00;\n"
"    float c0 = coeff * wz0;\n"
"    float c1 = coeff * wz1;\n"
"    uint base0 = evector_index(p, 0u, kk, ij);\n"
"    uint base1 = base0 + 3u;\n"
"    e1 += c0 * evector[base0] + c1 * evector[base1];\n"
"    e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];\n"
"    e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];\n"
"    ij++;\n"
"    coeff = wp0 * wt10;\n"
"    c0 = coeff * wz0;\n"
"    c1 = coeff * wz1;\n"
"    base0 = evector_index(p, 0u, kk, ij);\n"
"    base1 = base0 + 3u;\n"
"    e1 += c0 * evector[base0] + c1 * evector[base1];\n"
"    e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];\n"
"    e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];\n"
"\n"
"    ij = (uint)jtion1[o];\n"
"    const float wp1 = 1.0f - wp0;\n"
"    const float wt01 = 1.0f - wtion1[o];\n"
"    const float wt11 = 1.0f - wt01;\n"
"    coeff = wp1 * wt01;\n"
"    c0 = coeff * wz0;\n"
"    c1 = coeff * wz1;\n"
"    base0 = evector_index(p, 0u, kk, ij);\n"
"    base1 = base0 + 3u;\n"
"    e1 += c0 * evector[base0] + c1 * evector[base1];\n"
"    e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];\n"
"    e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];\n"
"    ij++;\n"
"    coeff = wp1 * wt11;\n"
"    c0 = coeff * wz0;\n"
"    c1 = coeff * wz1;\n"
"    base0 = evector_index(p, 0u, kk, ij);\n"
"    base1 = base0 + 3u;\n"
"    e1 += c0 * evector[base0] + c1 * evector[base1];\n"
"    e2 += c0 * evector[base0 + 1u] + c1 * evector[base1 + 1u];\n"
"    e3 += c0 * evector[base0 + 2u] + c1 * evector[base1 + 2u];\n"
"  }\n"
"  const float epsi = 0.25f * e1;\n"
"  const float etheta = 0.25f * e2;\n"
"  const float ezeta = 0.25f * e3;\n"
"\n"
"  device float *z = zion + m * (uint)GTC_NPARAM;\n"
"  device float *z0 = zion0 + m * (uint)GTC_NPARAM;\n"
"  if (p.irk == 1) {\n"
"    z0[GTC_Z_PSI] = z[GTC_Z_PSI];\n"
"    z0[GTC_Z_THETA] = z[GTC_Z_THETA];\n"
"    z0[GTC_Z_ZETA] = z[GTC_Z_ZETA];\n"
"    z0[GTC_Z_U] = z[GTC_Z_U];\n"
"    z0[GTC_Z_WEIGHT] = z[GTC_Z_WEIGHT];\n"
"  }\n"
"\n"
"  const float psi = z[GTC_Z_PSI];\n"
"  const float theta = z[GTC_Z_THETA];\n"
"  const float u = z[GTC_Z_U];\n"
"  const float weight = z[GTC_Z_WEIGHT];\n"
"  const float mu = z[GTC_Z_MU];\n"
"  const float r = sqrt(2.0f * psi);\n"
"  const float rinv = 1.0f / r;\n"
"  const int ii = clamp_int((int)floor((r - p.a0) * p.delr), 0, p.mpsi - 1);\n"
"  const float wp0 = (float)(ii + 1) - (r - p.a0) * p.delr;\n"
"  const float wp1 = 1.0f - wp0;\n"
"  const float tem = wp0 * temp_inv[ii] + wp1 * temp_inv[ii + 1];\n"
"  const float q = p.q0 + p.q1 * r * p.ainv + p.q2 * r * r * p.ainv * p.ainv;\n"
"  const float qinv = 1.0f / q;\n"
"  const float sint = sin(theta);\n"
"  const float cost = cos(theta);\n"
"  const float b = 1.0f / (1.0f + r * cost);\n"
"  const float g = 1.0f;\n"
"  const float gp = 0.0f;\n"
"  const float ri = 0.0f;\n"
"  const float rip = 0.0f;\n"
"  const float dbdp = -b * b * cost * rinv;\n"
"  const float dbdt = b * b * r * sint;\n"
"  const float dedb = p.cinv * (u * u * p.qion * b * p.cmratio + mu * mu);\n"
"  const float deni = 1.0f / (g * q + ri + u * (g * rip - ri * gp));\n"
"  const float upara = u * b * p.cmratio;\n"
"  const float energy = 0.5f * p.aion * upara * upara + mu * mu * b;\n"
"  float rfac = p.rw * (r - p.rc);\n"
"  rfac = rfac * rfac;\n"
"  rfac = rfac * rfac * rfac;\n"
"  rfac = exp(-rfac);\n"
"  const float kappa = ((energy * tem - 1.5f) * p.kappati + p.kappan) *\n"
"      (1.0f - p.sbound + p.sbound * rfac) * rinv;\n"
"\n"
"  float dptdp = epsi;\n"
"  float dptdt = etheta;\n"
"  float dptdz = ezeta - etheta * qinv;\n"
"  const float epara = -ezeta * b * q * deni;\n"
"  const float vdr = q * (ri * dptdz - g * dptdt) * deni;\n"
"  const float wdrive = vdr * kappa;\n"
"  const float wpara = epara * upara * p.qion * tem;\n"
"  const float wdrift = q * (g * dbdt * dptdp - g * dbdp * dptdt + ri * dbdp * dptdz) *\n"
"      deni * dedb * p.qion * tem;\n"
"  const float marker_weight = z0[GTC_Z_MU];\n"
"  const float wdot = (marker_weight - p.paranl * weight) * (wdrive + wpara + wdrift);\n"
"\n"
"  const float profile_vdr = p.has_vdrtmp != 0 ? vdrtmp[ii] : 0.0f;\n"
"  float pdot;\n"
"  float tdot;\n"
"  float zdot;\n"
"  float rdot;\n"
"  if (p.linear_orbit != 0) {\n"
"    pdot = q * (-dedb * dbdt) * deni;\n"
"    tdot = (upara * b + q * (dedb * dbdp)) * deni;\n"
"    zdot = upara * b * q * deni;\n"
"    rdot = -dedb * dbdt * deni;\n"
"  } else {\n"
"    dptdp = dptdp * p.nonlinear + p.gyroradius * p.flow0 * (1.0f - 0.5f * p.a * rinv);\n"
"    dptdt *= p.nonlinear;\n"
"    dptdz *= p.nonlinear;\n"
"    pdot = q * (-g * dedb * dbdt - g * dptdt + ri * dptdz) * deni - profile_vdr;\n"
"    tdot = (upara * b * (1.0f - q * gp * u) + q * g * (dedb * dbdp + dptdp)) * deni;\n"
"    zdot = (upara * b * q * (1.0f + rip * u) - q * ri * (dedb * dbdp + dptdp)) * deni;\n"
"    rdot = ((gp * u - 1.0f) * (dedb * dbdt + p.paranl * dptdt) -\n"
"        p.paranl * q * (1.0f + rip * u) * dptdz) * deni;\n"
"  }\n"
"\n"
"  z[GTC_Z_PSI] = max(1.0e-8f * p.psimax, z0[GTC_Z_PSI] + p.dtime * pdot);\n"
"  z[GTC_Z_THETA] = modulo_positive(z0[GTC_Z_THETA] + p.dtime * tdot, p.pi2);\n"
"  z[GTC_Z_ZETA] = modulo_positive(z0[GTC_Z_ZETA] + p.dtime * zdot, p.pi2);\n"
"  z[GTC_Z_U] = z0[GTC_Z_U] + p.dtime * rdot;\n"
"  z[GTC_Z_WEIGHT] = z0[GTC_Z_WEIGHT] + p.dtime * wdot;\n"
"\n"
"  if (p.write_wpi != 0) {\n"
"    wpi[0u * (uint)p.mimax + m] = vdr;\n"
"    wpi[1u * (uint)p.mimax + m] = energy;\n"
"    wpi[2u * (uint)p.mimax + m] = b;\n"
"  }\n"
"}\n";

static id<MTLDevice> gtc_metal_device;
static id<MTLCommandQueue> gtc_metal_queue;
static id<MTLLibrary> gtc_metal_library;
static id<MTLComputePipelineState> gtc_chargei_prepare_pipeline;
static id<MTLComputePipelineState> gtc_pushi_general_pipeline;
static int gtc_metal_initialized;
static int gtc_metal_enabled;
static int gtc_metal_reported_error;

static int gtc_env_truthy(const char *value) {
  if (!value || value[0] == '\0') return 0;
  if (strcmp(value, "0") == 0 || strcmp(value, "false") == 0 ||
      strcmp(value, "False") == 0 || strcmp(value, "off") == 0 ||
      strcmp(value, "OFF") == 0 || strcmp(value, "cpu") == 0 ||
      strcmp(value, "CPU") == 0) {
    return 0;
  }
  return 1;
}

static int gtc_env_requests_metal(void) {
  const char *gpu = getenv("GTC_GPU");
  if (gpu && gpu[0] != '\0') {
    return gtc_env_truthy(gpu);
  }
  return gtc_env_truthy(getenv("GTC_USE_GPU"));
}

static void gtc_metal_error_once(const char *message, NSError *error) {
  if (gtc_metal_reported_error) return;
  if (error) {
    fprintf(stderr, "gtc_c Metal backend disabled: %s: %s\n", message,
            [[error localizedDescription] UTF8String]);
  } else {
    fprintf(stderr, "gtc_c Metal backend disabled: %s\n", message);
  }
  gtc_metal_reported_error = 1;
}

static int gtc_metal_init(void) {
  if (gtc_metal_initialized) return gtc_metal_enabled;
  gtc_metal_initialized = 1;
  gtc_metal_enabled = 0;

  if (!gtc_env_requests_metal()) return 0;

  @autoreleasepool {
    gtc_metal_device = MTLCreateSystemDefaultDevice();
    if (!gtc_metal_device) {
      gtc_metal_error_once("no default Metal device", nil);
      return 0;
    }
    gtc_metal_queue = [gtc_metal_device newCommandQueue];
    if (!gtc_metal_queue) {
      gtc_metal_error_once("could not create command queue", nil);
      return 0;
    }

    NSError *error = nil;
    NSString *source = [NSString stringWithUTF8String:gtc_metal_source];
    MTLCompileOptions *options = [MTLCompileOptions new];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    options.fastMathEnabled = NO;
#pragma clang diagnostic pop
    gtc_metal_library = [gtc_metal_device newLibraryWithSource:source
                                                       options:options
                                                         error:&error];
    if (!gtc_metal_library) {
      gtc_metal_error_once("could not compile runtime kernels", error);
      return 0;
    }
  }

  gtc_metal_enabled = 1;
  return 1;
}

int gtc_gpu_enabled(void) {
  return gtc_metal_init();
}

const char *gtc_gpu_backend_name(void) {
  return gtc_metal_init() ? "metal" : "disabled";
}

static id<MTLComputePipelineState> gtc_metal_pipeline(const char *name) {
  if (!gtc_metal_init()) return nil;

  id<MTLComputePipelineState> __strong *slot = NULL;
  if (strcmp(name, "chargei_prepare_kernel") == 0) {
    slot = &gtc_chargei_prepare_pipeline;
  } else if (strcmp(name, "pushi_general_kernel") == 0) {
    slot = &gtc_pushi_general_pipeline;
  } else {
    gtc_metal_error_once("unknown kernel name", nil);
    return nil;
  }
  if (*slot) return *slot;

  @autoreleasepool {
    NSString *function_name = [NSString stringWithUTF8String:name];
    id<MTLFunction> function = [gtc_metal_library newFunctionWithName:function_name];
    if (!function) {
      gtc_metal_error_once("could not find kernel function", nil);
      return nil;
    }
    NSError *error = nil;
    id<MTLComputePipelineState> pipeline =
        [gtc_metal_device newComputePipelineStateWithFunction:function error:&error];
    if (!pipeline) {
      gtc_metal_error_once("could not create compute pipeline", error);
      return nil;
    }
    *slot = pipeline;
  }
  return *slot;
}

static id<MTLBuffer> gtc_metal_buffer_copy(const void *data, size_t length) {
  if (length == 0) length = 1;
  return [gtc_metal_device newBufferWithBytes:data
                                       length:length
                                      options:MTLResourceStorageModeShared];
}

static id<MTLBuffer> gtc_metal_buffer_empty(size_t length) {
  if (length == 0) length = 1;
  return [gtc_metal_device newBufferWithLength:length
                                       options:MTLResourceStorageModeShared];
}

static id<MTLBuffer> gtc_metal_buffer_host(void *data, size_t length,
                                           int initialize_fallback, int *direct) {
  if (length == 0) length = 1;
  if (direct) *direct = 0;
  if (!data) return gtc_metal_buffer_empty(length);

  id<MTLBuffer> buffer =
      [gtc_metal_device newBufferWithBytesNoCopy:data
                                          length:length
                                         options:MTLResourceStorageModeShared
                                     deallocator:nil];
  if (buffer) {
    if (direct) *direct = 1;
    return buffer;
  }
  return initialize_fallback ? gtc_metal_buffer_copy(data, length)
                             : gtc_metal_buffer_empty(length);
}

static void gtc_metal_copy_back_if_needed(id<MTLBuffer> buffer, int direct,
                                          void *host, size_t length) {
  if (!direct && buffer && host && length > 0) {
    memcpy(host, [buffer contents], length);
  }
}

static int gtc_metal_dispatch(id<MTLComputePipelineState> pipeline,
                              id<MTLBuffer> __strong *buffers, int buffer_count,
                              NSUInteger work_items) {
  if (work_items == 0) return 1;
  id<MTLCommandBuffer> command_buffer = [gtc_metal_queue commandBuffer];
  id<MTLComputeCommandEncoder> encoder = [command_buffer computeCommandEncoder];
  if (!command_buffer || !encoder) {
    gtc_metal_error_once("could not create command buffer", nil);
    return 0;
  }

  [encoder setComputePipelineState:pipeline];
  for (int i = 0; i < buffer_count; i++) {
    if (!buffers[i]) {
      gtc_metal_error_once("could not allocate Metal buffer", nil);
      return 0;
    }
    [encoder setBuffer:buffers[i] offset:0 atIndex:(NSUInteger)i];
  }

  NSUInteger width = [pipeline maxTotalThreadsPerThreadgroup];
  if (width > 256) width = 256;
  if (width < 1) width = 1;
  MTLSize grid_size = MTLSizeMake(work_items, 1, 1);
  MTLSize group_size = MTLSizeMake(width, 1, 1);
  [encoder dispatchThreads:grid_size threadsPerThreadgroup:group_size];
  [encoder endEncoding];
  [command_buffer commit];
  [command_buffer waitUntilCompleted];

  if ([command_buffer status] == MTLCommandBufferStatusError) {
    gtc_metal_error_once("command buffer failed", [command_buffer error]);
    return 0;
  }
  return 1;
}

int gtc_gpu_chargei_prepare(GtcState *s, const GtcReal *delt, GtcReal delr,
                            GtcReal delz, GtcReal smu_inv, GtcReal pi2_inv,
                            int single_zeta_cell) {
  if (!s || s->mi <= 0 || !gtc_metal_init()) return 0;

  @autoreleasepool {
    id<MTLComputePipelineState> pipeline = gtc_metal_pipeline("chargei_prepare_kernel");
    if (!pipeline) return 0;

    const size_t particle_count = (size_t)s->mi;
    const size_t zion_bytes = particle_count * (size_t)GTC_NPARAM * sizeof(GtcReal);
    const size_t larmor_particle_bytes = particle_count * 4u * sizeof(GtcReal);
    const size_t larmor_particle_int_bytes = particle_count * 4u * sizeof(int);
    const size_t radial_bytes = (size_t)(s->p.mpsi + 1) * sizeof(GtcReal);
    const size_t radial_int_bytes = (size_t)(s->p.mpsi + 1) * sizeof(int);
    const size_t gyro_bytes = (size_t)s->mgrid * 4u * sizeof(GtcReal);
    const size_t particle_int_bytes = particle_count * sizeof(int);
    const size_t particle_real_bytes = particle_count * sizeof(GtcReal);

    GtcMetalChargeParams params = {
      .mi = s->mi,
      .mpsi = s->p.mpsi,
      .mzeta = s->mzeta,
      .single_zeta_cell = single_zeta_cell,
      .delr = delr,
      .delz = delz,
      .smu_inv = smu_inv,
      .pi2_inv = pi2_inv,
      .zetamin = s->zetamin,
      .a0 = s->p.a0,
      .a1 = s->p.a1
    };

    int direct[15] = {0};
    id<MTLBuffer> __strong buffers[15];
    buffers[0] = gtc_metal_buffer_host(s->zion, zion_bytes, 1, &direct[0]);
    buffers[1] = gtc_metal_buffer_host(s->pgyro, gyro_bytes, 1, &direct[1]);
    buffers[2] = gtc_metal_buffer_host(s->tgyro, gyro_bytes, 1, &direct[2]);
    buffers[3] = gtc_metal_buffer_host(s->mtheta, radial_int_bytes, 1, &direct[3]);
    buffers[4] = gtc_metal_buffer_host(s->igrid, radial_int_bytes, 1, &direct[4]);
    buffers[5] = gtc_metal_buffer_host(s->qtinv, radial_bytes, 1, &direct[5]);
    buffers[6] = gtc_metal_buffer_host((void *)delt, radial_bytes, 1, &direct[6]);
    buffers[7] = gtc_metal_buffer_host(s->kzion, particle_int_bytes, 0, &direct[7]);
    buffers[8] = gtc_metal_buffer_host(s->wzion, particle_real_bytes, 0, &direct[8]);
    buffers[9] = gtc_metal_buffer_host(s->wpion, larmor_particle_bytes, 0, &direct[9]);
    buffers[10] = gtc_metal_buffer_host(s->jtion0, larmor_particle_int_bytes, 0, &direct[10]);
    buffers[11] = gtc_metal_buffer_host(s->jtion1, larmor_particle_int_bytes, 0, &direct[11]);
    buffers[12] = gtc_metal_buffer_host(s->wtion0, larmor_particle_bytes, 0, &direct[12]);
    buffers[13] = gtc_metal_buffer_host(s->wtion1, larmor_particle_bytes, 0, &direct[13]);
    buffers[14] = gtc_metal_buffer_copy(&params, sizeof(params));

    if (!gtc_metal_dispatch(pipeline, buffers, 15, (NSUInteger)particle_count)) return 0;

    gtc_metal_copy_back_if_needed(buffers[7], direct[7], s->kzion, particle_int_bytes);
    gtc_metal_copy_back_if_needed(buffers[8], direct[8], s->wzion, particle_real_bytes);
    gtc_metal_copy_back_if_needed(buffers[9], direct[9], s->wpion, larmor_particle_bytes);
    gtc_metal_copy_back_if_needed(buffers[10], direct[10], s->jtion0, larmor_particle_int_bytes);
    gtc_metal_copy_back_if_needed(buffers[11], direct[11], s->jtion1, larmor_particle_int_bytes);
    gtc_metal_copy_back_if_needed(buffers[12], direct[12], s->wtion0, larmor_particle_bytes);
    gtc_metal_copy_back_if_needed(buffers[13], direct[13], s->wtion1, larmor_particle_bytes);
  }
  return 1;
}

static int gtc_gpu_pushi_general_impl(GtcState *s, const GtcReal *temp_inv,
                                      const GtcReal *vdrtmp, int has_vdrtmp,
                                      int linear_orbit, int copy_wpi_back,
                                      GtcReal delr, GtcReal pi2,
                                      GtcReal psimax, GtcReal cmratio,
                                      GtcReal cinv, GtcReal vthi,
                                      GtcReal ainv, GtcReal sbound,
                                      GtcReal dtime) {
  if (!s || s->mi <= 0 || !gtc_metal_init()) return 0;

  @autoreleasepool {
    id<MTLComputePipelineState> pipeline = gtc_metal_pipeline("pushi_general_kernel");
    if (!pipeline) return 0;

    const size_t particle_count = (size_t)s->mi;
    const size_t active_particle_bytes = particle_count * (size_t)GTC_NPARAM * sizeof(GtcReal);
    const size_t wpi_bytes = 3u * (size_t)s->mimax * sizeof(GtcReal);
    const size_t field_bytes =
        (size_t)(s->mzeta + 1) * (size_t)s->mgrid * 3u * sizeof(GtcReal);
    const size_t particle_int_bytes = particle_count * sizeof(int);
    const size_t particle_real_bytes = particle_count * sizeof(GtcReal);
    const size_t larmor_particle_bytes = particle_count * 4u * sizeof(GtcReal);
    const size_t larmor_particle_int_bytes = particle_count * 4u * sizeof(int);
    const size_t radial_bytes = (size_t)(s->p.mpsi + 1) * sizeof(GtcReal);
    const GtcReal dummy_vdrtmp = 0.0f;

    GtcMetalPushiParams params = {
      .mi = s->mi,
      .mimax = s->mimax,
      .mpsi = s->p.mpsi,
      .mzeta = s->mzeta,
      .irk = s->irk,
      .linear_orbit = linear_orbit,
      .has_vdrtmp = has_vdrtmp,
      .write_wpi = copy_wpi_back,
      .delr = delr,
      .pi2 = pi2,
      .psimax = psimax,
      .cmratio = cmratio,
      .cinv = cinv,
      .vthi = vthi,
      .ainv = ainv,
      .sbound = sbound,
      .dtime = dtime,
      .a = s->p.a,
      .a0 = s->p.a0,
      .q0 = s->p.q0,
      .q1 = s->p.q1,
      .q2 = s->p.q2,
      .qion = s->p.qion,
      .aion = s->p.aion,
      .kappati = s->p.kappati,
      .kappan = s->p.kappan,
      .paranl = s->p.paranl,
      .nonlinear = s->p.nonlinear,
      .flow0 = s->p.flow0,
      .rc = s->p.rc,
      .rw = s->p.rw,
      .gyroradius = s->gyroradius
    };

    int direct[14] = {0};
    id<MTLBuffer> __strong buffers[14];
    buffers[0] = gtc_metal_buffer_host(s->zion, active_particle_bytes, 1, &direct[0]);
    buffers[1] = gtc_metal_buffer_host(s->zion0, active_particle_bytes, 1, &direct[1]);
    buffers[2] = copy_wpi_back ? gtc_metal_buffer_host(s->wpi, wpi_bytes, 0, &direct[2])
                               : gtc_metal_buffer_empty(sizeof(GtcReal));
    buffers[3] = gtc_metal_buffer_host(s->evector, field_bytes, 1, &direct[3]);
    buffers[4] = gtc_metal_buffer_host(s->kzion, particle_int_bytes, 1, &direct[4]);
    buffers[5] = gtc_metal_buffer_host(s->jtion0, larmor_particle_int_bytes, 1, &direct[5]);
    buffers[6] = gtc_metal_buffer_host(s->jtion1, larmor_particle_int_bytes, 1, &direct[6]);
    buffers[7] = gtc_metal_buffer_host(s->wzion, particle_real_bytes, 1, &direct[7]);
    buffers[8] = gtc_metal_buffer_host(s->wpion, larmor_particle_bytes, 1, &direct[8]);
    buffers[9] = gtc_metal_buffer_host(s->wtion0, larmor_particle_bytes, 1, &direct[9]);
    buffers[10] = gtc_metal_buffer_host(s->wtion1, larmor_particle_bytes, 1, &direct[10]);
    buffers[11] = gtc_metal_buffer_host((void *)temp_inv, radial_bytes, 1, &direct[11]);
    buffers[12] = has_vdrtmp ?
        gtc_metal_buffer_host((void *)vdrtmp, radial_bytes, 1, &direct[12]) :
        gtc_metal_buffer_copy(&dummy_vdrtmp, sizeof(dummy_vdrtmp));
    buffers[13] = gtc_metal_buffer_copy(&params, sizeof(params));

    if (!gtc_metal_dispatch(pipeline, buffers, 14, (NSUInteger)particle_count)) return 0;

    gtc_metal_copy_back_if_needed(buffers[0], direct[0], s->zion, active_particle_bytes);
    gtc_metal_copy_back_if_needed(buffers[1], direct[1], s->zion0, active_particle_bytes);
    if (copy_wpi_back) gtc_metal_copy_back_if_needed(buffers[2], direct[2], s->wpi, wpi_bytes);
  }
  return 1;
}

int gtc_gpu_pushi_general(GtcState *s, const GtcReal *temp_inv,
                          const GtcReal *vdrtmp, int has_vdrtmp,
                          int linear_orbit, GtcReal delr, GtcReal pi2,
                          GtcReal psimax, GtcReal cmratio, GtcReal cinv,
                          GtcReal vthi, GtcReal ainv, GtcReal sbound,
                          GtcReal dtime) {
  return gtc_gpu_pushi_general_impl(s, temp_inv, vdrtmp, has_vdrtmp, linear_orbit, 1,
                                   delr, pi2, psimax, cmratio, cinv, vthi, ainv,
                                   sbound, dtime);
}

int gtc_gpu_pushi_linear_orbit(GtcState *s, const GtcReal *temp_inv,
                               GtcReal delr, GtcReal pi2, GtcReal psimax,
                               GtcReal cmratio, GtcReal cinv, GtcReal vthi,
                               GtcReal ainv, GtcReal sbound, GtcReal dtime) {
  const GtcReal dummy_vdrtmp = 0.0f;
  return gtc_gpu_pushi_general_impl(s, temp_inv, &dummy_vdrtmp, 0, 1, 0,
                                   delr, pi2, psimax, cmratio, cinv, vthi, ainv,
                                   sbound, dtime);
}
