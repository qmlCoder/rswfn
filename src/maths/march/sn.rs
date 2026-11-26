//! Surface Nets（Dual Contouring 的简化版）等值面提取。
//!
//! 每个与等值面相交的立方体生成 **一个** 顶点（取该立方体 12 条棱上交点位置的平均），
//! 对每条「符号发生变化的网格棱」生成一个四边形（拆成两个三角形）连接它周围的 4 个立方体。
//! 顶点天然按立方体共享，输出水密、流形，且三角形数量比 Marching Cubes / Tetrahedra 都少。

use super::{DIDX, EDGES, grid_index};

/// 棱 `(pa, pb)` 上等值点（`va`、`vb` 必不相等）
#[inline]
fn crossing(pa: [f64; 3], pb: [f64; 3], va: f64, vb: f64, isov: f64) -> [f64; 3] {
    let t = (isov - va) / (vb - va);
    [
        pa[0] + t * (pb[0] - pa[0]),
        pa[1] + t * (pb[1] - pa[1]),
        pa[2] + t * (pb[2] - pa[2]),
    ]
}

/// 对单个等值面 `isov` 提取三角形（扁平顶点表，每 3 个一组）。
fn extract_iso(shape: &[usize; 3], grids: &[[f64; 3]], values: &[f64], isov: f64) -> Vec<[f64; 3]> {
    let [nx, ny, nz] = *shape;
    let mut tris: Vec<[f64; 3]> = Vec::new();
    if nx < 2 || ny < 2 || nz < 2 {
        return tris;
    }
    let (cx, cy, cz) = (nx - 1, ny - 1, nz - 1);
    // 立方体 (a, b, c) 在 cell_v 中的下标
    let cidx = |a: usize, b: usize, c: usize| (a * cy + b) * cz + c;

    // 第一遍：给每个相交立方体算一个顶点
    let mut cell_v = vec![usize::MAX; cx * cy * cz];
    let mut verts: Vec<[f64; 3]> = Vec::new();
    for i in 0..cx {
        for j in 0..cy {
            for k in 0..cz {
                let mut p = [[0.0f64; 3]; 8];
                let mut v = [0.0f64; 8];
                let mut nan = false;
                for c in 0..8 {
                    let [di, dj, dk] = DIDX[c];
                    let gi = grid_index(i + di, j + dj, k + dk, ny, nz);
                    p[c] = grids[gi];
                    v[c] = values[gi];
                    nan |= v[c].is_nan();
                }
                if nan {
                    continue;
                }
                let vmin = v.iter().cloned().fold(f64::INFINITY, f64::min);
                let vmax = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                if vmin > isov || vmax < isov {
                    continue;
                }
                // 12 条棱上交点位置的平均
                let mut sum = [0.0f64; 3];
                let mut cnt = 0.0;
                for e in EDGES.iter() {
                    let a = e[0] as usize;
                    let b = e[1] as usize;
                    if (v[a] > isov) != (v[b] > isov) {
                        let q = crossing(p[a], p[b], v[a], v[b], isov);
                        sum[0] += q[0];
                        sum[1] += q[1];
                        sum[2] += q[2];
                        cnt += 1.0;
                    }
                }
                if cnt == 0.0 {
                    continue;
                }
                verts.push([sum[0] / cnt, sum[1] / cnt, sum[2] / cnt]);
                cell_v[cidx(i, j, k)] = verts.len() - 1;
            }
        }
    }

    // 第二遍：每条符号变化的网格棱 -> 一个四边形。
    // emit 把四边形 (a,b,c,d) 拆成两个三角形追加到输出（法线朝「值减小」的一侧）。
    let emit = |a: usize, b: usize, c: usize, d: usize, out: &mut Vec<[f64; 3]>| {
        if a == usize::MAX || b == usize::MAX || c == usize::MAX || d == usize::MAX {
            return; // 相邻立方体里有 NaN 被跳过
        }
        out.push(verts[a]);
        out.push(verts[b]);
        out.push(verts[c]);
        out.push(verts[a]);
        out.push(verts[c]);
        out.push(verts[d]);
    };

    // x 方向的网格棱
    for i in 0..cx {
        for j in 1..cy {
            for k in 1..cz {
                let s0 = values[grid_index(i, j, k, ny, nz)];
                let s1 = values[grid_index(i + 1, j, k, ny, nz)];
                if s0.is_nan() || s1.is_nan() || (s0 > isov) == (s1 > isov) {
                    continue;
                }
                let c00 = cell_v[cidx(i, j - 1, k - 1)];
                let c10 = cell_v[cidx(i, j, k - 1)];
                let c11 = cell_v[cidx(i, j, k)];
                let c01 = cell_v[cidx(i, j - 1, k)];
                if s0 > isov {
                    emit(c00, c10, c11, c01, &mut tris);
                } else {
                    emit(c00, c01, c11, c10, &mut tris);
                }
            }
        }
    }

    // y 方向的网格棱
    for j in 0..cy {
        for i in 1..cx {
            for k in 1..cz {
                let s0 = values[grid_index(i, j, k, ny, nz)];
                let s1 = values[grid_index(i, j + 1, k, ny, nz)];
                if s0.is_nan() || s1.is_nan() || (s0 > isov) == (s1 > isov) {
                    continue;
                }
                let a = cell_v[cidx(i - 1, j, k - 1)];
                let b = cell_v[cidx(i - 1, j, k)];
                let c = cell_v[cidx(i, j, k)];
                let d = cell_v[cidx(i, j, k - 1)];
                if s0 > isov {
                    emit(a, b, c, d, &mut tris);
                } else {
                    emit(a, d, c, b, &mut tris);
                }
            }
        }
    }

    // z 方向的网格棱
    for k in 0..cz {
        for i in 1..cx {
            for j in 1..cy {
                let s0 = values[grid_index(i, j, k, ny, nz)];
                let s1 = values[grid_index(i, j, k + 1, ny, nz)];
                if s0.is_nan() || s1.is_nan() || (s0 > isov) == (s1 > isov) {
                    continue;
                }
                let a = cell_v[cidx(i - 1, j - 1, k)];
                let b = cell_v[cidx(i, j - 1, k)];
                let c = cell_v[cidx(i, j, k)];
                let d = cell_v[cidx(i - 1, j, k)];
                if s0 > isov {
                    emit(a, b, c, d, &mut tris);
                } else {
                    emit(a, d, c, b, &mut tris);
                }
            }
        }
    }

    tris
}

/// Surface Nets 提取等值面（接口与 [`super::mc::march_cube_v1`] 一致）。
///
/// 返回 `(顶点表, 类型表)`：类型 `+1` 表示 `isov` 层，`-1` 表示 `-isov` 层。
/// 注意：为与 [`super::mc`] / [`super::mt`] 保持一致的「三角形 soup」输出，
/// 这里把四边形展开成了扁平顶点；若需要带索引的网格可自行按 cell 顺序复用顶点。
pub fn surface_nets_v1(
    shape: &[usize; 3],
    grids: &Vec<[f64; 3]>,
    values: &Vec<f64>,
    isov: f64,
) -> (Vec<[f64; 3]>, Vec<isize>) {
    let mut pverts = extract_iso(shape, grids, values, isov);
    let mut nverts = extract_iso(shape, grids, values, -isov);
    nverts.reverse(); // 负层反转绕向

    let mut types = vec![1isize; pverts.len()];
    types.extend(vec![-1isize; nverts.len()]);
    pverts.extend(nverts);
    (pverts, types)
}
