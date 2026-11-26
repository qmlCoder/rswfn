//! Marching Tetrahedra（行进四面体）等值面提取。
//!
//! 把每个立方体按 Kuhn 分解切成 6 个四面体（都共用主对角线 `0-6`）。
//! 一个四面体与等值面的交只有「1/3」或「2/2」两种，没有歧义，因此生成的网格天然水密。
//! 代价是三角形数量比 Marching Cubes 多、光滑面略糙。

use super::{Cube, DIDX, grid_index};

/// 立方体的 6 个四面体（Kuhn 分解，共用主对角线 `0-6`，用立方体顶点号 0..8）
const TETS: [[usize; 4]; 6] = [
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
];

/// 一个四面体的 6 条棱（用四面体局部顶点号 0..4）
const TET_EDGES: [[usize; 2]; 6] = [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]];

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// 处理一个四面体，把三角形的顶点按「每 3 个一组」追加到 `out`。
///
/// 所有三角形的法线统一朝「`val < isov`」一侧，保证同一等值面内部绕向一致。
fn tet_triangles(p: &[[f64; 3]; 4], v: &[f64; 4], isov: f64, out: &mut Vec<[f64; 3]>) {
    let inside = [v[0] > isov, v[1] > isov, v[2] > isov, v[3] > isov];
    let nin = inside.iter().filter(|b| **b).count();
    if nin == 0 || nin == 4 {
        return;
    }

    // 内部 / 外部顶点质心，用来确定「朝外」方向
    let mut in_c = [0.0; 3];
    let mut out_c = [0.0; 3];
    let mut in_n = 0.0;
    let mut out_n = 0.0;
    for i in 0..4 {
        if inside[i] {
            in_c[0] += p[i][0];
            in_c[1] += p[i][1];
            in_c[2] += p[i][2];
            in_n += 1.0;
        } else {
            out_c[0] += p[i][0];
            out_c[1] += p[i][1];
            out_c[2] += p[i][2];
            out_n += 1.0;
        }
    }
    let dir = [
        out_c[0] / out_n - in_c[0] / in_n,
        out_c[1] / out_n - in_c[1] / in_n,
        out_c[2] / out_n - in_c[2] / in_n,
    ];

    // 四面体棱 (a, b) 上的交点
    let point = |a: usize, b: usize| -> [f64; 3] {
        let t = (isov - v[a]) / (v[b] - v[a]);
        [
            p[a][0] + t * (p[b][0] - p[a][0]),
            p[a][1] + t * (p[b][1] - p[a][1]),
            p[a][2] + t * (p[b][2] - p[a][2]),
        ]
    };

    // 按朝外方向压入一个三角形
    let mut push = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
        let n = cross(sub(b, a), sub(c, a));
        if dot(n, dir) < 0.0 {
            out.push(a);
            out.push(c);
            out.push(b);
        } else {
            out.push(a);
            out.push(b);
            out.push(c);
        }
    };

    // 收集交叉棱
    let cross_edges: Vec<(usize, usize)> = TET_EDGES
        .iter()
        .filter(|e| inside[e[0]] != inside[e[1]])
        .map(|e| (e[0], e[1]))
        .collect();

    if cross_edges.len() == 3 {
        let p0 = point(cross_edges[0].0, cross_edges[0].1);
        let p1 = point(cross_edges[1].0, cross_edges[1].1);
        let p2 = point(cross_edges[2].0, cross_edges[2].1);
        push(p0, p1, p2);
    } else {
        // 4 条交叉棱 -> 四边形（拆成两个三角形）。
        // 相邻的交叉棱共享一个四面体顶点，顺着共享顶点把 4 条棱排成环。
        let n = cross_edges.len();
        let mut order = vec![0usize];
        let mut used = vec![false; n];
        used[0] = true;
        while order.len() < n {
            let last = cross_edges[*order.last().unwrap()];
            let mut found = false;
            for k in 0..n {
                if used[k] {
                    continue;
                }
                let e = cross_edges[k];
                let share = last.0 == e.0 || last.0 == e.1 || last.1 == e.0 || last.1 == e.1;
                if share {
                    order.push(k);
                    used[k] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }
        let mut quad = Vec::with_capacity(4);
        for &idx in &order {
            let e = cross_edges[idx];
            quad.push(point(e.0, e.1));
        }
        if quad.len() == 4 {
            push(quad[0], quad[1], quad[2]);
            push(quad[0], quad[2], quad[3]);
        }
    }
}

/// 对单个等值面 `isov` 提取三角形（扁平顶点表，每 3 个一组）。
fn extract_iso(shape: &[usize; 3], grids: &[[f64; 3]], values: &[f64], isov: f64) -> Vec<[f64; 3]> {
    let [nx, ny, nz] = *shape;
    let mut out: Vec<[f64; 3]> = Vec::new();
    if nx < 2 || ny < 2 || nz < 2 {
        return out;
    }

    for i in 0..nx - 1 {
        for j in 0..ny - 1 {
            for k in 0..nz - 1 {
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

                for tet in TETS.iter() {
                    let tp = [p[tet[0]], p[tet[1]], p[tet[2]], p[tet[3]]];
                    let tv = [v[tet[0]], v[tet[1]], v[tet[2]], v[tet[3]]];
                    tet_triangles(&tp, &tv, isov, &mut out);
                }
            }
        }
    }
    out
}

/// Marching Tetrahedra 提取等值面（接口与 [`super::mc::march_cube_v1`] 一致）。
///
/// 返回 `(顶点表, 类型表)`：类型 `+1` 表示 `isov` 层，`-1` 表示 `-isov` 层。
pub fn march_tetra_v1(
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

/// 在给定区域内用 Marching Tetrahedra 提取等值面，并做平滑处理（`func` 为标量场）。
///
/// MT 的分割方向固定，直接提取出来的表面会显得皱；这里在提取后把网格焊接起来做
/// Taubin 平滑，再重新投影回真实等值面，在保持形状的同时让结果变平滑
/// （代价：多几次 `func` 求值）。
pub fn march_tetra_v2<F>(
    func: F,
    cube: Cube,
    isov: f64,
    step: f64,
    xyzs: Option<&Vec<[f64; 3]>>,
    limit: Option<f64>,
) -> (Vec<[f64; 3]>, Vec<isize>)
where
    F: Fn(&[f64; 3]) -> f64 + Sync + Send,
{
    let (verts, types) = super::extract_region(&func, cube, isov, step, xyzs, limit, march_tetra_v1);
    // 焊接 + Taubin 平滑 + 重新投影回等值面
    let verts = super::smooth_isosurface(&verts, &types, isov, &func, 6);
    (verts, types)
}
