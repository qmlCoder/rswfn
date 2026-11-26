//! 终端里的分子结构预览。
//!
//! 分子坐标先按当前视角旋转，再做正交投影，光栅化成"每个字符格两个像素"的
//! 半块图像：`▀` 的上半块用前景色、下半块用背景色，于是一个字符格能显示
//! 两种颜色，横向 1 像素、纵向 2 像素。球体按漫反射 + 高光着色，化学键按
//! 圆柱着色，远处原子向背景雾化（景深），所以纯文本终端也能看出立体感。
//!
//! 同一份 [`Frame`] 既给 TUI 生成 [`Line`]，也给 CLI 生成 ANSI 序列，
//! 两者视觉完全一致。颜色来自 `datas::elems` 里的 CPK 配色。

use std::io::IsTerminal;
use std::path::Path;

use ratatui::style::{Color, Style as TuiStyle};
use ratatui::text::{Line, Span};

use crate::base::Mole;
use crate::datas::consts::BOHR;
use crate::datas::elems;

/// 视图背景色。抗锯齿和景深雾化都向它混合，TUI 里同时作为面板底色，
/// 这样终端自带背景色不会透出来影响球体边缘。
pub const BG: Rgb = Rgb(16, 16, 22);

/// 真彩色 RGB。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rgb(pub u8, pub u8, pub u8);

impl Rgb {
    /// 转为 ratatui 颜色
    /// 按终端色彩能力映射成 ratatui 颜色
    pub fn color_in(self, depth: Depth) -> Color {
        match depth {
            Depth::True => Color::Rgb(self.0, self.1, self.2),
            Depth::Ansi256 => Color::Indexed(self.ansi256()),
            Depth::Mono => Color::Reset,
        }
    }

    pub fn color(self) -> Color {
        self.color_in(Depth::True)
    }

    /// 最接近的 xterm-256 调色板索引（6×6×6 色立方 + 24 级灰阶）
    pub fn ansi256(self) -> u8 {
        const STEPS: [u16; 6] = [0, 95, 135, 175, 215, 255];

        fn nearest(value: u8) -> usize {
            let mut best = 0;
            let mut best_d = u16::MAX;
            for (i, step) in STEPS.iter().enumerate() {
                let d = (value as u16).abs_diff(*step);
                if d < best_d {
                    best_d = d;
                    best = i;
                }
            }
            best
        }

        fn dist2(a: Rgb, b: Rgb) -> u32 {
            let d = |x: u8, y: u8| (x as i32 - y as i32).unsigned_abs().pow(2);
            d(a.0, b.0) + d(a.1, b.1) + d(a.2, b.2)
        }

        let ri = nearest(self.0);
        let gi = nearest(self.1);
        let bi = nearest(self.2);
        let cube = Rgb(STEPS[ri] as u8, STEPS[gi] as u8, STEPS[bi] as u8);
        let cube_idx = 16 + 36 * ri as u8 + 6 * gi as u8 + bi as u8;

        // 灰阶 232..=255 对应 8, 18, …, 238
        let mean = (self.0 as u32 + self.1 as u32 + self.2 as u32) / 3;
        let gray_idx = if mean > 238 { 23 } else { (mean.saturating_sub(3) / 10).min(23) as u8 };
        let step = 8 + 10 * gray_idx;
        let gray = Rgb(step, step, step);

        if dist2(self, gray) < dist2(self, cube) { 232 + gray_idx } else { cube_idx }
    }

    /// 解析 `#RRGGBB`
    pub fn from_hex(hex: &str) -> Option<Rgb> {
        let hex = hex.trim().trim_start_matches('#');
        if hex.len() != 6 {
            return None;
        }
        let r = u8::from_str_radix(&hex[0..2], 16).ok()?;
        let g = u8::from_str_radix(&hex[2..4], 16).ok()?;
        let b = u8::from_str_radix(&hex[4..6], 16).ok()?;
        Some(Rgb(r, g, b))
    }

    /// `self` 向 `other` 混合，`t=0` 返回 self，`t=1` 返回 other
    fn mix(self, other: Rgb, t: f64) -> Rgb {
        let t = t.clamp(0.0, 1.0);
        let lerp = |a: u8, b: u8| (a as f64 + (b as f64 - a as f64) * t).round().clamp(0.0, 255.0) as u8;
        Rgb(lerp(self.0, other.0), lerp(self.1, other.1), lerp(self.2, other.2))
    }

    /// 整体亮度缩放
    fn scale(self, k: f64) -> Rgb {
        let f = |v: u8| (v as f64 * k).round().clamp(0.0, 255.0) as u8;
        Rgb(f(self.0), f(self.1), f(self.2))
    }

    /// 叠加白光（高光）
    fn add(self, k: f64) -> Rgb {
        let f = |v: u8| (v as f64 + k * 255.0).round().clamp(0.0, 255.0) as u8;
        Rgb(f(self.0), f(self.1), f(self.2))
    }

    /// 感知亮度 0..1
    fn lum(self) -> f64 {
        (0.2126 * self.0 as f64 + 0.7152 * self.1 as f64 + 0.0722 * self.2 as f64) / 255.0
    }
}

/// 元素符号对应的 CPK 颜色，未知元素退化为灰色。
fn element_color(nuc: usize) -> Rgb {
    elems::ELEMS.iter().find(|elem| elem.idx == nuc as isize).and_then(|elem| Rgb::from_hex(&elem.color)).unwrap_or(Rgb(180, 180, 180))
}

/// 用于绘制的分子快照：只保留坐标、半径、颜色和成键关系。
///
/// `Mole` 里的 `Atoms` 带缓存与锁，不便跨线程克隆；预览只读一份扁平数据。
#[derive(Debug, Clone)]
pub struct Molecule {
    /// 元素符号
    pub syms: Vec<String>,
    /// 原子坐标（Bohr，保持文件里的原始朝向）
    pub coords: Vec<[f64; 3]>,
    /// 共价半径（Bohr）
    pub rads: Vec<f64>,
    /// CPK 颜色
    pub colors: Vec<Rgb>,
    /// 成键关系（0 基原子序号对）
    pub bonds: Vec<[usize; 2]>,
    /// 主轴方向：三行分别是屏幕 x、y、深度的世界坐标方向
    pub axes: [[f64; 3]; 3],
}

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

impl Default for Molecule {
    fn default() -> Self {
        Self {
            syms: Vec::new(),
            coords: Vec::new(),
            rads: Vec::new(),
            colors: Vec::new(),
            bonds: Vec::new(),
            axes: IDENTITY,
        }
    }
}

impl Molecule {
    /// 从 `Mole` 提取预览所需的扁平常量数据
    pub fn from_mole(mole: &Mole) -> Self {
        let atoms = mole.atoms();
        let coords = atoms.xyzs().clone();
        let center = centroid_of(&coords);
        let axes = principal_axes(&coords, center);
        Self {
            syms: atoms.syms(),
            coords,
            rads: atoms.rads(),
            colors: atoms.nucs().iter().map(|&nuc| element_color(nuc)).collect(),
            bonds: atoms.bonds(),
            axes,
        }
    }

    /// 原子数
    pub fn len(&self) -> usize {
        self.coords.len()
    }

    /// 是否没有原子
    pub fn is_empty(&self) -> bool {
        self.coords.is_empty()
    }

    /// 几何中心
    pub fn centroid(&self) -> [f64; 3] {
        centroid_of(&self.coords)
    }
}

fn centroid_of(coords: &[[f64; 3]]) -> [f64; 3] {
    let n = coords.len().max(1) as f64;
    let mut c = [0.0; 3];
    for xyz in coords {
        for k in 0..3 {
            c[k] += xyz[k];
        }
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

/// 主轴方向（坐标协方差矩阵的特征向量），按摊开程度从大到小排列。
///
/// 分子文件里的朝向是任意的，直接投影经常看到"侧着"的苯环。主成分分析把
/// 摊得最开的两个方向摆到屏幕 x / y，最窄的方向留给深度，平面分子就正对镜头。
fn principal_axes(coords: &[[f64; 3]], c: [f64; 3]) -> [[f64; 3]; 3] {
    if coords.len() < 3 {
        return IDENTITY;
    }
    let mut a = [[0.0f64; 3]; 3];
    for p in coords {
        let v = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += v[i] * v[j];
            }
        }
    }

    // Jacobi 旋转，把 3×3 对称矩阵对角化（不用引入线代依赖）
    let mut v = IDENTITY;
    for _ in 0..32 {
        let (mut p, mut q, mut best) = (0usize, 1usize, a[0][1].abs());
        if a[0][2].abs() > best {
            p = 0;
            q = 2;
            best = a[0][2].abs();
        }
        if a[1][2].abs() > best {
            p = 1;
            q = 2;
            best = a[1][2].abs();
        }
        let scale = (a[0][0].abs() + a[1][1].abs() + a[2][2].abs()).max(1e-30);
        if best <= 1e-12 * scale {
            break;
        }
        let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
        let cos = 1.0 / (1.0 + t * t).sqrt();
        let sin = t * cos;
        for k in 0..3 {
            if k == p || k == q {
                continue;
            }
            let akp = a[k][p];
            let akq = a[k][q];
            let np = cos * akp - sin * akq;
            let nq = sin * akp + cos * akq;
            a[k][p] = np;
            a[p][k] = np;
            a[k][q] = nq;
            a[q][k] = nq;
        }
        let (app, aqq, apq) = (a[p][p], a[q][q], a[p][q]);
        a[p][p] = cos * cos * app - 2.0 * sin * cos * apq + sin * sin * aqq;
        a[q][q] = sin * sin * app + 2.0 * sin * cos * apq + cos * cos * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for k in 0..3 {
            let vkp = v[k][p];
            let vkq = v[k][q];
            v[k][p] = cos * vkp - sin * vkq;
            v[k][q] = sin * vkp + cos * vkq;
        }
    }

    // 特征值从大到小排序，对应的特征向量就是三个主轴
    let mut order = [0usize, 1, 2];
    order.sort_by(|&i, &j| a[j][j].partial_cmp(&a[i][i]).unwrap_or(std::cmp::Ordering::Equal));
    let mut axes = [[v[0][order[0]], v[1][order[0]], v[2][order[0]]], [v[0][order[1]], v[1][order[1]], v[2][order[1]]], [0.0; 3]];

    // 特征向量的正负号是任意的，用三次矩把朝向定下来，保证同一分子每次一样
    for axis in axes.iter_mut().take(2) {
        let mut skew = 0.0;
        for p in coords {
            let d = (p[0] - c[0]) * axis[0] + (p[1] - c[1]) * axis[1] + (p[2] - c[2]) * axis[2];
            skew += d * d * d;
        }
        if skew < 0.0 {
            for k in 0..3 {
                axis[k] = -axis[k];
            }
        }
    }
    // 右手系：深度方向 = x × y
    axes[2] = [
        axes[0][1] * axes[1][2] - axes[0][2] * axes[1][1],
        axes[0][2] * axes[1][0] - axes[0][0] * axes[1][2],
        axes[0][0] * axes[1][1] - axes[0][1] * axes[1][0],
    ];
    axes
}

/// 绘制模型。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Style {
    /// 球棍模型
    #[default]
    BallStick,
    /// 空间填充（CPK）模型
    SpaceFilling,
    /// 线框模型
    Wireframe,
}

impl Style {
    /// 中文名（用于界面标题）
    pub fn name(self) -> &'static str {
        match self {
            Style::BallStick => "球棍",
            Style::SpaceFilling => "空间填充",
            Style::Wireframe => "线框",
        }
    }

    /// 循环切换
    pub fn next(self) -> Style {
        match self {
            Style::BallStick => Style::SpaceFilling,
            Style::SpaceFilling => Style::Wireframe,
            Style::Wireframe => Style::BallStick,
        }
    }

    /// 解析命令行取值: ball | space | wire
    pub fn parse(value: &str) -> Option<Style> {
        match value {
            "ball" | "ball-and-stick" | "sticks" => Some(Style::BallStick),
            "space" | "spacefill" | "cpk" => Some(Style::SpaceFilling),
            "wire" | "wireframe" | "line" => Some(Style::Wireframe),
            _ => None,
        }
    }
}

/// 视角与绘制参数。
#[derive(Debug, Clone)]
pub struct View {
    /// 视角旋转（四元数 w,x,y,z）。按屏幕轴累加，不做俯仰角限制，
    /// 行为类似 three.js 的 TrackballControls。
    rot: [f64; 4],
    /// 缩放倍数（1.0 为自动铺满）
    pub zoom: f64,
    /// 绘制模型
    pub style: Style,
    /// 是否在原子球上写元素符号
    pub labels: bool,
    /// 是否启用景深雾化
    pub depth_cue: bool,
}

fn quat_axis(axis: [f64; 3], angle: f64) -> [f64; 4] {
    let (sin, cos) = (angle * 0.5).sin_cos();
    [cos, axis[0] * sin, axis[1] * sin, axis[2] * sin]
}

fn quat_mul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    let [aw, ax, ay, az] = a;
    let [bw, bx, by, bz] = b;
    [
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ]
}

fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if norm < 1e-12 { [1.0, 0.0, 0.0, 0.0] } else { [q[0] / norm, q[1] / norm, q[2] / norm, q[3] / norm] }
}

fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [w, x, y, z] = q;
    // t = 2 * (qv × v)，v' = v + w*t + qv × t
    let tx = 2.0 * (y * v[2] - z * v[1]);
    let ty = 2.0 * (z * v[0] - x * v[2]);
    let tz = 2.0 * (x * v[1] - y * v[0]);
    [v[0] + w * tx + (y * tz - z * ty), v[1] + w * ty + (z * tx - x * tz), v[2] + w * tz + (x * ty - y * tx)]
}

impl Default for View {
    fn default() -> Self {
        // 与旧默认视角 yaw=0.35、pitch=0.21 对应（先绕 y、后绕 x）
        let rot = quat_mul(quat_axis([1.0, 0.0, 0.0], 0.21), quat_axis([0.0, 1.0, 0.0], 0.35));
        Self {
            rot: quat_normalize(rot),
            zoom: 1.0,
            style: Style::BallStick,
            labels: true,
            depth_cue: true,
        }
    }
}

impl View {
    /// Trackball 旋转（弧度增量）：绕当前屏幕的 x/y 轴累加，
    /// 越过极点后仍可继续旋转，不会因欧拉角限制而卡住。
    pub fn rotate(&mut self, dx: f64, dy: f64) {
        let delta = quat_mul(quat_axis([1.0, 0.0, 0.0], dy), quat_axis([0.0, 1.0, 0.0], dx));
        self.rot = quat_normalize(quat_mul(delta, self.rot));
    }

    /// 直接设置绝对视角（弧度）；不限制 pitch，允许翻转
    pub fn angles(&mut self, yaw: f64, pitch: f64) {
        let rot = quat_mul(quat_axis([1.0, 0.0, 0.0], pitch), quat_axis([0.0, 1.0, 0.0], yaw));
        self.rot = quat_normalize(rot);
    }

    /// 缩放（倍数，>1 放大）
    pub fn zoom_by(&mut self, factor: f64) {
        self.zoom = (self.zoom * factor).clamp(0.2, 8.0);
    }

    /// 切换绘制模型
    pub fn next_style(&mut self) {
        self.style = self.style.next();
    }

    /// 每个原子的绘制半径（Bohr），自适应缩放要把球本身也算进外廓。
    fn ball_weights(&self, mol: &Molecule) -> Vec<f64> {
        (0..mol.len())
            .map(|i| {
                let ang = mol.rads[i] / BOHR;
                let radius_ang = match self.style {
                    Style::BallStick => 0.09 + 0.26 * ang,
                    Style::SpaceFilling => 2.0 * ang,
                    // 线框没有球，也要留一点余量给线宽和抗锯齿
                    Style::Wireframe => 0.06,
                };
                radius_ang * BOHR
            })
            .collect()
    }

    /// 旋转不变的包围球半径（质心到原子 + 绘制半径）。
    ///
    /// 不要用投影后的二维包围盒：分子一转，投影宽度就变，
    /// 自适应缩放下分子会看起来忽大忽小，侧视时还会被放得过大。
    fn bounding_radius(&self, mol: &Molecule) -> f64 {
        let c = mol.centroid();
        let ball_w = self.ball_weights(mol);
        let mut radius: f64 = 1e-9;
        for i in 0..mol.len() {
            let v = [mol.coords[i][0] - c[0], mol.coords[i][1] - c[1], mol.coords[i][2] - c[2]];
            let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            radius = radius.max(len + ball_w[i]);
        }
        radius
    }

    /// 自适应缩放：只由包围球、画布尺寸和 zoom 决定，与当前视角无关。
    fn fit_scale(&self, mol: &Molecule, width: u16, height: u16) -> f64 {
        let w = width.max(1) as f64;
        let h = height.max(1) as f64;
        0.44 * w.min(2.0 * h) / self.bounding_radius(mol) * self.zoom
    }

    /// 渲染为字符帧
    ///
    /// `width` / `height` 是字符格数，实际像素网格为 `width × (2 * height)`。
    pub fn render(&self, mol: &Molecule, width: u16, height: u16) -> Frame {
        let width = width.max(1);
        let height = height.max(1);
        let w = width as usize;
        let h = height as usize;
        let mut frame = Frame {
            width,
            height,
            cells: vec![Cell { top: BG, bottom: BG, over: None }; w * h],
        };
        if mol.is_empty() {
            return frame;
        }

        // ---- 投影：转到主轴坐标系 → 四元数 Trackball 旋转 ----
        let c = mol.centroid();
        let n = mol.len();
        let ball_w = self.ball_weights(mol);
        let mut pts = vec![[0.0f64; 3]; n];
        let mut zmin = f64::INFINITY;
        let mut zmax = f64::NEG_INFINITY;
        for i in 0..n {
            let v = [mol.coords[i][0] - c[0], mol.coords[i][1] - c[1], mol.coords[i][2] - c[2]];
            let ax = dot(mol.axes[0], v);
            let ay = dot(mol.axes[1], v);
            let az = dot(mol.axes[2], v);
            let [x1, y1, z2] = quat_rotate(self.rot, [ax, ay, az]);
            pts[i] = [x1, y1, z2];
            zmin = zmin.min(z2);
            zmax = zmax.max(z2);
        }

        // 自适应缩放：用旋转不变的包围球，旋转时大小保持稳定
        let scale = self.fit_scale(mol, width, height);
        let cx = w as f64 / 2.0;
        let cy = h as f64;
        let screen: Vec<[f64; 2]> = pts.iter().map(|p| [cx + p[0] * scale, cy - p[1] * scale]).collect();

        // ---- 景深：把 z 映射成 0.62..1.0 的可见度 ----
        let depth_span = (zmax - zmin).max(1e-9);
        let fog = |z: f64| -> f64 { if !self.depth_cue { 1.0 } else { 0.62 + 0.38 * ((z - zmin) / depth_span).clamp(0.0, 1.0) } };

        // ---- 半径（像素） ----
        let ball_px = |i: usize| (ball_w[i] * scale).max(if self.style == Style::Wireframe { 0.9 } else { 1.1 });
        let bond_px = match self.style {
            Style::BallStick => (0.14 * BOHR * scale).max(0.9),
            Style::Wireframe => (0.07 * BOHR * scale).max(0.9),
            Style::SpaceFilling => 0.0,
        };

        // ---- 画笔排序：远处的先画（画家算法） ----
        #[derive(Clone, Copy)]
        enum Item {
            Atom(usize),
            Bond(usize),
        }
        let mut items: Vec<(f64, Item)> = Vec::with_capacity(n + mol.bonds.len());
        for i in 0..n {
            items.push((pts[i][2], Item::Atom(i)));
        }
        if bond_px > 0.0 {
            for (k, [i, j]) in mol.bonds.iter().enumerate() {
                items.push(((pts[*i][2] + pts[*j][2]) / 2.0, Item::Bond(k)));
            }
        }
        items.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let mut pix = Pixmap::new(w, 2 * h);
        for (_, item) in items {
            match item {
                Item::Atom(i) => {
                    let color = mol.colors.get(i).copied().unwrap_or(Rgb(180, 180, 180));
                    pix.sphere(screen[i][0], screen[i][1], ball_px(i), color, fog(pts[i][2]));
                }
                Item::Bond(k) => {
                    let [i, j] = mol.bonds[k];
                    let ca = mol.colors.get(i).copied().unwrap_or(Rgb(180, 180, 180));
                    let cb = mol.colors.get(j).copied().unwrap_or(Rgb(180, 180, 180));
                    // 键画到球面为止，否则圆柱会盖住球心，看起来像"糖葫芦"
                    let (a, b) = (screen[i], screen[j]);
                    let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
                    let len = (dx * dx + dy * dy).sqrt().max(1e-9);
                    let trim_a = (ball_px(i) * 0.85).min(len * 0.45);
                    let trim_b = (ball_px(j) * 0.85).min(len * 0.45);
                    let a2 = [a[0] + dx / len * trim_a, a[1] + dy / len * trim_a];
                    let b2 = [b[0] - dx / len * trim_b, b[1] - dy / len * trim_b];
                    pix.bond(a2, b2, bond_px, ca, cb, fog(pts[i][2]), fog(pts[j][2]));
                }
            }
        }

        // ---- 元素符号：由远及近覆盖，近处原子优先 ----
        if self.labels {
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|&i, &j| pts[i][2].partial_cmp(&pts[j][2]).unwrap_or(std::cmp::Ordering::Equal));
            for i in order {
                let radius = ball_px(i);
                if radius < 2.4 {
                    continue;
                }
                let chars: Vec<char> = mol.syms[i].chars().collect();
                if chars.is_empty() || chars.len() > 3 {
                    continue;
                }
                let row = (screen[i][1] / 2.0).floor();
                if row < 0.0 || row >= height as f64 {
                    continue;
                }
                let start = (screen[i][0] - chars.len() as f64 / 2.0).round();
                let row = row as usize;
                for (k, ch) in chars.iter().enumerate() {
                    let col = start + k as f64;
                    if col < 0.0 || col >= width as f64 {
                        continue;
                    }
                    let col = col as usize;
                    let bg = pix.get(col, 2 * row);
                    let fg = if bg.lum() > 0.52 { Rgb(12, 12, 16) } else { Rgb(246, 246, 250) };
                    frame.cells[row * w + col].over = Some((*ch, fg));
                }
            }
        }

        for row in 0..h {
            for col in 0..w {
                let cell = &mut frame.cells[row * w + col];
                cell.top = pix.get(col, 2 * row);
                cell.bottom = pix.get(col, 2 * row + 1);
            }
        }
        frame
    }
}

/// 终端色彩能力
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Depth {
    /// 24 位真彩色（`38;2;r;g;b`）
    #[default]
    True,
    /// 256 色（`38;5;n`）
    Ansi256,
    /// 无色，用 ASCII 灰度字符
    Mono,
}

impl Depth {
    /// 解析显式指定的色彩深度
    pub fn parse(value: &str) -> Option<Depth> {
        match value.trim().to_ascii_lowercase().as_str() {
            "truecolor" | "true" | "24bit" | "rgb" => Some(Depth::True),
            "256" | "256color" | "ansi256" | "ansi" => Some(Depth::Ansi256),
            "mono" | "none" | "ascii" | "gray" => Some(Depth::Mono),
            _ => None,
        }
    }

    /// 根据环境变量推断终端能承受的色彩深度。
    ///
    /// 关键点：真彩色必须由终端自己支持。不支持时若仍输出 `38;2;r;g;b`，
    /// 这些参数会被当成普通 SGR / 调色板索引重新解释，整幅图就变成随机色带
    /// （Xshell 这类只有 256 色的终端就是这样）。
    pub fn detect() -> Depth {
        fn var(name: &str) -> String {
            std::env::var(name).unwrap_or_default().to_ascii_lowercase()
        }

        let colorterm = var("COLORTERM");
        if colorterm.contains("truecolor") || colorterm.contains("24bit") {
            return Depth::True;
        }

        let term = var("TERM");
        if term.contains("truecolor") || term.contains("direct") || term.contains("24bit") {
            return Depth::True;
        }

        // 已知支持 24 位色的终端，但它们不一定在 TERM/COLORTERM 里声明
        let program = var("TERM_PROGRAM");
        if matches!(program.as_str(), "iterm.app" | "wezterm" | "vscode" | "kitty" | "hyper" | "tabby" | "ghostty" | "warpterminal") {
            return Depth::True;
        }
        // Windows Terminal 与新版 conhost 都支持 24 位色，而它们的 TERM 可能为空
        if !var("WT_SESSION").is_empty() {
            return Depth::True;
        }

        if term == "dumb" {
            return Depth::Mono;
        }
        if term.is_empty() {
            return if cfg!(windows) { Depth::True } else { Depth::Mono };
        }
        Depth::Ansi256
    }
}

/// 单个字符格：上半块颜色、下半块颜色、可选覆盖字符（元素符号）。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Cell {
    top: Rgb,
    bottom: Rgb,
    over: Option<(char, Rgb)>,
}

impl Cell {
    /// 渲染成 (字符, 前景色, 背景色)
    fn paint(&self) -> (char, Rgb, Rgb) {
        match self.over {
            Some((ch, fg)) => (ch, fg, self.top),
            None if self.top == self.bottom => ('█', self.top, self.top),
            None => ('▀', self.top, self.bottom),
        }
    }
}

/// 无色模式下按亮度映射的 ASCII 灰度字符
fn mono_char(cell: &Cell) -> char {
    const RAMP: [char; 10] = [' ', '.', ':', '-', '=', '+', '*', '#', '%', '@'];
    if let Some((ch, _)) = cell.over {
        return ch;
    }
    // 背景色（很暗）一律留白，否则整幅图会铺满 '.' 看不出分子轮廓
    let lum = cell.top.mix(cell.bottom, 0.5).lum();
    let idx = if lum < 0.10 { 0 } else { (((lum - 0.10) / 0.90 * 9.0).round() as usize).clamp(1, 9) };
    RAMP[idx]
}

/// 渲染好的字符帧。
#[derive(Debug, Clone)]
pub struct Frame {
    /// 字符列数
    pub width: u16,
    /// 字符行数
    pub height: u16,
    cells: Vec<Cell>,
}

impl Frame {
    /// 生成 ratatui 文本行（TUI 用），按终端色彩能力着色
    pub fn lines_in(&self, depth: Depth) -> Vec<Line<'static>> {
        if depth == Depth::Mono {
            return (0..self.height as usize)
                .map(|row| {
                    let mut text = String::with_capacity(self.width as usize);
                    for col in 0..self.width as usize {
                        text.push(mono_char(&self.cells[row * self.width as usize + col]));
                    }
                    Line::from(text)
                })
                .collect();
        }
        let w = self.width as usize;
        (0..self.height as usize)
            .map(|row| {
                let mut spans: Vec<Span<'static>> = Vec::new();
                let mut text = String::new();
                let mut style: Option<(Rgb, Rgb)> = None;
                for col in 0..w {
                    let (ch, fg, bg) = self.cells[row * w + col].paint();
                    if style != Some((fg, bg)) {
                        if !text.is_empty() {
                            let (f, b) = style.expect("样式存在");
                            spans.push(Span::styled(std::mem::take(&mut text), TuiStyle::default().fg(f.color_in(depth)).bg(b.color_in(depth))));
                        }
                        style = Some((fg, bg));
                    }
                    text.push(ch);
                }
                if !text.is_empty() {
                    let (f, b) = style.expect("样式存在");
                    spans.push(Span::styled(text, TuiStyle::default().fg(f.color_in(depth)).bg(b.color_in(depth))));
                }
                Line::from(spans)
            })
            .collect()
    }

    /// 真彩色版本（测试与默认使用）
    pub fn lines(&self) -> Vec<Line<'static>> {
        self.lines_in(Depth::True)
    }

    /// 生成带 ANSI 转义序列的文本（CLI 用），按终端色彩能力着色
    pub fn ansi_in(&self, depth: Depth) -> String {
        let w = self.width as usize;
        let mut out = String::new();
        for row in 0..self.height as usize {
            let mut style: Option<(Rgb, Rgb)> = None;
            for col in 0..w {
                let (ch, fg, bg) = self.cells[row * w + col].paint();
                if style != Some((fg, bg)) {
                    match depth {
                        Depth::True => out.push_str(&format!("\x1b[38;2;{};{};{};48;2;{};{};{}m", fg.0, fg.1, fg.2, bg.0, bg.1, bg.2)),
                        _ => out.push_str(&format!("\x1b[38;5;{};48;5;{}m", fg.ansi256(), bg.ansi256())),
                    }
                    style = Some((fg, bg));
                }
                out.push(ch);
            }
            out.push_str("\x1b[0m\n");
        }
        out
    }

    /// 真彩色版本（测试与默认使用）
    pub fn ansi(&self) -> String {
        self.ansi_in(Depth::True)
    }

    /// 无色输出：按亮度映射成 ASCII 灰度字符（重定向到文件或无色终端时用）
    pub fn plain(&self) -> String {
        let w = self.width as usize;
        let mut out = String::new();
        for row in 0..self.height as usize {
            for col in 0..w {
                out.push(mono_char(&self.cells[row * w + col]));
            }
            out.push('\n');
        }
        out
    }
}

// ===== 光栅化 =====

const LIGHT: [f64; 3] = [-0.45, 0.55, 0.70];

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-9);
    [v[0] / len, v[1] / len, v[2] / len]
}

/// 漫反射 + 高光着色
fn shade(n: [f64; 3], base: Rgb) -> Rgb {
    let l = normalize(LIGHT);
    let h = normalize([l[0], l[1], l[2] + 1.0]);
    let diff = (n[0] * l[0] + n[1] * l[1] + n[2] * l[2]).max(0.0);
    let spec = (n[0] * h[0] + n[1] * h[1] + n[2] * h[2]).max(0.0).powf(32.0) * 0.18;
    base.scale(0.30 + 0.70 * diff).add(spec)
}

/// 像素画布（不透明，底色为 [`BG`]）
struct Pixmap {
    w: usize,
    h: usize,
    px: Vec<Rgb>,
}

impl Pixmap {
    fn new(w: usize, h: usize) -> Self {
        Self { w, h, px: vec![BG; w * h] }
    }

    fn get(&self, x: usize, y: usize) -> Rgb {
        self.px[y.min(self.h - 1) * self.w + x.min(self.w - 1)]
    }

    fn blend(&mut self, x: i64, y: i64, color: Rgb, alpha: f64) {
        if x < 0 || y < 0 || x >= self.w as i64 || y >= self.h as i64 {
            return;
        }
        let idx = y as usize * self.w + x as usize;
        self.px[idx] = self.px[idx].mix(color, alpha);
    }

    /// 着色球体
    fn sphere(&mut self, cx: f64, cy: f64, r: f64, base: Rgb, fog: f64) {
        let r = r.max(0.6);
        let x0 = (cx - r - 1.0).floor() as i64;
        let x1 = (cx + r + 1.0).ceil() as i64;
        let y0 = (cy - r - 1.0).floor() as i64;
        let y1 = (cy + r + 1.0).ceil() as i64;
        for y in y0..=y1 {
            for x in x0..=x1 {
                let dx = x as f64 + 0.5 - cx;
                let dy = y as f64 + 0.5 - cy;
                let d2 = dx * dx + dy * dy;
                let d = d2.sqrt();
                let cover = (r + 0.5 - d).clamp(0.0, 1.0);
                if cover <= 0.0 {
                    continue;
                }
                let nz = (1.0 - d2 / (r * r)).max(0.0).sqrt();
                let color = shade([dx / r, -dy / r, nz], base).mix(BG, 1.0 - fog);
                self.blend(x, y, color, cover);
            }
        }
    }

    /// 着色化学键（圆柱），`ca`→`cb` 沿键方向渐变
    fn bond(&mut self, a: [f64; 2], b: [f64; 2], r: f64, ca: Rgb, cb: Rgb, fa: f64, fb: f64) {
        let r = r.max(0.35);
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let len2 = dx * dx + dy * dy;
        let x0 = (a[0].min(b[0]) - r - 1.0).floor() as i64;
        let x1 = (a[0].max(b[0]) + r + 1.0).ceil() as i64;
        let y0 = (a[1].min(b[1]) - r - 1.0).floor() as i64;
        let y1 = (a[1].max(b[1]) + r + 1.0).ceil() as i64;
        for y in y0..=y1 {
            for x in x0..=x1 {
                let px = x as f64 + 0.5;
                let py = y as f64 + 0.5;
                let t = if len2 > 1e-9 { (((px - a[0]) * dx + (py - a[1]) * dy) / len2).clamp(0.0, 1.0) } else { 0.0 };
                let ox = px - (a[0] + dx * t);
                let oy = py - (a[1] + dy * t);
                let d = (ox * ox + oy * oy).sqrt();
                let cover = (r + 0.5 - d).clamp(0.0, 1.0);
                if cover <= 0.0 {
                    continue;
                }
                let side = (d / r).clamp(0.0, 1.0);
                let (nx, ny) = if d > 1e-9 { (ox / d * side, oy / d * side) } else { (0.0, 0.0) };
                let nz = (1.0 - side * side).max(0.0).sqrt();
                let color = shade([nx, -ny, nz], ca.mix(cb, t)).mix(BG, 1.0 - (fa + (fb - fa) * t));
                // 细线在网格上采样时容易断成虚线，抬一点覆盖率让笔画连续
                self.blend(x, y, color, cover.powf(0.7));
            }
        }
    }
}

// ===== CLI 支持 =====

/// 颜色输出模式
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ColorMode {
    /// 输出到终端时上色，重定向时退化为灰度字符
    #[default]
    Auto,
    /// 总是输出 ANSI 真彩色
    Always,
    /// 从不输出 ANSI
    Never,
}

impl ColorMode {
    /// 解析命令行取值
    pub fn parse(value: &str) -> Option<ColorMode> {
        match value {
            "auto" => Some(ColorMode::Auto),
            "always" => Some(ColorMode::Always),
            "never" => Some(ColorMode::Never),
            _ => None,
        }
    }

    fn enabled(self) -> bool {
        match self {
            ColorMode::Auto => std::io::stdout().is_terminal(),
            ColorMode::Always => true,
            ColorMode::Never => false,
        }
    }
}

/// `rswfn view` 的绘制选项
#[derive(Debug, Clone)]
pub struct RenderOptions {
    /// 画布尺寸（列 × 行），None 表示按终端大小自动选择
    pub size: Option<(u16, u16)>,
    /// 绘制模型
    pub style: Style,
    /// 初始视角（度）：绕 y 轴、绕 x 轴
    pub rotate: (f64, f64),
    /// 缩放倍数
    pub zoom: f64,
    /// 是否显示元素符号
    pub labels: bool,
    /// 是否启用景深
    pub depth_cue: bool,
    /// 颜色模式
    pub color: ColorMode,
    /// 输出 JSON（原子坐标与成键）而不是字符画
    pub json: bool,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            size: None,
            style: Style::BallStick,
            rotate: (20.0, 12.0),
            zoom: 1.0,
            labels: true,
            depth_cue: true,
            color: ColorMode::Auto,
            json: false,
        }
    }
}

impl RenderOptions {
    /// 画布尺寸：显式指定优先，否则用终端大小（留一行给提示），再不行 80×32
    fn canvas(&self) -> (u16, u16) {
        if let Some(size) = self.size {
            return (size.0.max(4), size.1.max(2));
        }
        match ratatui::crossterm::terminal::size() {
            Ok((cols, rows)) => (cols.max(4), rows.saturating_sub(1).max(2)),
            Err(_) => (80, 32),
        }
    }
}

/// 读取文件并构造预览用分子。
pub fn load(input: &Path) -> Result<Molecule, String> {
    let mole = crate::service::load_mole(input).or_else(|_| crate::service::load_any_mole(input))?;
    Ok(Molecule::from_mole(&mole))
}

/// 把分子导出为 JSON（1 基序号，坐标单位 Å）。
pub fn to_json(mol: &Molecule, input: &Path) -> String {
    let atoms: Vec<serde_json::Value> = (0..mol.len())
        .map(|i| {
            serde_json::json!({
                "index": i + 1,
                "symbol": mol.syms.get(i).cloned().unwrap_or_default(),
                "x": mol.coords[i][0] / BOHR,
                "y": mol.coords[i][1] / BOHR,
                "z": mol.coords[i][2] / BOHR,
                "radius": mol.rads[i] / BOHR,
            })
        })
        .collect();
    let bonds: Vec<[usize; 2]> = mol.bonds.iter().map(|[i, j]| [i + 1, j + 1]).collect();
    serde_json::to_string_pretty(&serde_json::json!({
        "file": input.display().to_string(),
        "natoms": mol.len(),
        "atoms": atoms,
        "bonds": bonds,
    }))
    .expect("预览 JSON 序列化失败")
}

/// `rswfn view` 的实现
pub fn run(input: &Path, opts: &RenderOptions) -> Result<(), String> {
    let mol = load(input)?;
    if mol.is_empty() {
        return Err(format!("{} 中没有原子坐标", input.display()));
    }
    if opts.json {
        println!("{}", to_json(&mol, input));
        return Ok(());
    }
    let mut view = View {
        zoom: opts.zoom,
        style: opts.style,
        labels: opts.labels,
        depth_cue: opts.depth_cue,
        ..View::default()
    };
    view.angles(opts.rotate.0.to_radians(), opts.rotate.1.to_radians());
    let (width, height) = opts.canvas();
    let frame = view.render(&mol, width, height);
    // 按终端能力选择色彩深度：不支持真彩色的终端（如只有 256 色的 Xshell）
    // 若硬输出 38;2;… 会被错解成随机色带，必须降级到 256 色或纯 ASCII。
    let text = if !opts.color.enabled() {
        frame.plain()
    } else {
        match Depth::detect() {
            Depth::Mono => frame.plain(),
            depth => frame.ansi_in(depth),
        }
    };
    print!("{text}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 水分子（Å → Bohr），半径取近似共价半径
    fn water() -> Molecule {
        let b = BOHR;
        Molecule {
            syms: vec!["O".to_string(), "H".to_string(), "H".to_string()],
            coords: vec![[0.0, 0.0, 0.0], [0.96 * b, 0.0, 0.0], [-0.24 * b, 0.93 * b, 0.0]],
            rads: vec![0.63 * b, 0.32 * b, 0.32 * b],
            colors: vec![Rgb(229, 0, 0), Rgb(204, 204, 204), Rgb(204, 204, 204)],
            bonds: vec![[0, 1], [0, 2]],
            axes: IDENTITY,
        }
    }

    /// 统计非空白字符（灰度输出里背景是空格）
    fn ink(frame: &Frame) -> usize {
        frame.plain().chars().filter(|ch| !ch.is_whitespace()).count()
    }

    #[test]
    fn renders_something() {
        let frame = View::default().render(&water(), 60, 20);
        assert_eq!(frame.width, 60);
        assert_eq!(frame.height, 20);
        assert!(ink(&frame) > 100, "画布上应该有可见的分子");
    }

    #[test]
    fn empty_molecule_is_blank() {
        let frame = View::default().render(&Molecule::default(), 40, 10);
        assert_eq!(ink(&frame), 0);
    }

    #[test]
    fn labels_show_element_symbols() {
        let view = View { labels: true, ..View::default() };
        let text: String = view.render(&water(), 80, 30).lines().iter().map(|line| line.to_string()).collect();
        assert!(text.contains('O'), "应显示氧的元素符号");
        assert!(text.contains('H'), "应显示氢的元素符号");
    }

    #[test]
    fn labels_can_be_disabled() {
        let view = View { labels: false, ..View::default() };
        let text: String = view.render(&water(), 80, 30).lines().iter().map(|line| line.to_string()).collect();
        assert!(!text.contains('O'));
    }

    #[test]
    fn ansi_uses_truecolor() {
        let ansi = View::default().render(&water(), 40, 12).ansi();
        // 前景色和背景色写在同一条 SGR 里：38;2;r;g;b;48;2;r;g;b
        assert!(ansi.contains("\x1b[38;2;"));
        assert!(ansi.contains(";48;2;"));
        assert!(ansi.ends_with("\x1b[0m\n"));
    }

    #[test]
    fn zoom_changes_ink() {
        let base = View::default().render(&water(), 60, 20);
        let zoomed = View { zoom: 2.0, ..View::default() }.render(&water(), 60, 20);
        assert!(ink(&zoomed) > ink(&base), "放大后像素应更多");
    }

    #[test]
    fn styles_differ() {
        let ball = View::default().render(&water(), 60, 20).plain();
        let wire = View {
            style: Style::Wireframe,
            ..View::default()
        }
        .render(&water(), 60, 20)
        .plain();
        let space = View {
            style: Style::SpaceFilling,
            ..View::default()
        }
        .render(&water(), 60, 20)
        .plain();
        assert_ne!(ball, wire);
        assert_ne!(ball, space);
        assert_ne!(wire, space);
    }

    #[test]
    fn rotation_changes_image() {
        let a = View::default().render(&water(), 60, 20).plain();
        let mut view = View::default();
        view.angles(1.9, -0.7);
        let b = view.render(&water(), 60, 20).plain();
        assert_ne!(a, b);
    }

    #[test]
    fn trackball_rotation_passes_poles() {
        let mut view = View::default();
        view.angles(0.0, 1.4); // 旧实现会在这里附近被夹住
        let a = view.render(&water(), 60, 20).plain();
        view.rotate(0.0, 0.4); // 越过 90°
        let b = view.render(&water(), 60, 20).plain();
        assert_ne!(a, b, "越过 90° 后仍应继续旋转");
        view.rotate(0.0, 0.4);
        let c = view.render(&water(), 60, 20).plain();
        assert_ne!(b, c, "越过极点后不应卡住");

        // 连续上下翻转多圈仍能渲染
        for _ in 0..100 {
            view.rotate(0.0, 0.5);
        }
        assert!(ink(&view.render(&water(), 60, 20)) > 0, "连续翻转后仍应可见");

        // 缩放仍然保留上下限
        view.zoom_by(1000.0);
        assert!(view.zoom <= 8.0);
        view.zoom_by(1e-6);
        assert!(view.zoom >= 0.2);
    }

    #[test]
    fn auto_scale_is_rotation_invariant() {
        let mol = water();
        let mut a = View::default();
        a.angles(0.0, 0.0);
        let mut b = View::default();
        b.angles(1.3, -0.4);
        let sa = a.fit_scale(&mol, 80, 30);
        let sb = b.fit_scale(&mol, 80, 30);
        assert!((sa - sb).abs() < 1e-12, "自适应缩放不应随视角变化（否则旋转时忽大忽小）");

        let zoomed = View { zoom: 2.0, ..View::default() };
        assert!((zoomed.fit_scale(&mol, 80, 30) - 2.0 * sa).abs() < 1e-12, "zoom 应线性影响缩放");
    }

    #[test]
    fn quat_rotation_preserves_length() {
        let q = quat_normalize(quat_mul(quat_axis([1.0, 0.0, 0.0], 1.3), quat_axis([0.0, 1.0, 0.0], -0.7)));
        let v = [0.3, -1.2, 2.5];
        let r = quat_rotate(q, v);
        let len = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        assert!((len(r) - len(v)).abs() < 1e-9);
    }

    #[test]
    fn parses_style_and_color_mode() {
        assert_eq!(Style::parse("space"), Some(Style::SpaceFilling));
        assert_eq!(Style::parse("wire"), Some(Style::Wireframe));
        assert_eq!(Style::parse("nope"), None);
        assert_eq!(ColorMode::parse("always"), Some(ColorMode::Always));
        assert_eq!(ColorMode::parse("nope"), None);
    }

    #[test]
    fn hex_colors_parse() {
        assert_eq!(Rgb::from_hex("#8E8E8E"), Some(Rgb(142, 142, 142)));
        assert_eq!(Rgb::from_hex("00ff00"), Some(Rgb(0, 255, 0)));
        assert_eq!(Rgb::from_hex("#fff"), None);
    }

    #[test]
    fn ansi256_quantizes_to_palette() {
        assert_eq!(Rgb(0, 0, 0).ansi256(), 16);
        assert_eq!(Rgb(255, 255, 255).ansi256(), 231);
        assert_eq!(Rgb(255, 0, 0).ansi256(), 196);
        // 灰色应当落在灰阶段（232..=255），而不是 6×6×6 色立方
        let gray = Rgb(128, 128, 128).ansi256();
        assert!((232..=255).contains(&gray), "灰色应落在灰阶段，实际 {gray}");
    }

    #[test]
    fn content_stays_inside_canvas() {
        let b = BOHR;
        // 苯环 + 氢，坐标取自常见几何，足够撑满画布
        let mut syms = Vec::new();
        let mut coords = Vec::new();
        let mut rads = Vec::new();
        let mut bonds = Vec::new();
        for k in 0..6 {
            let a = std::f64::consts::TAU * k as f64 / 6.0;
            syms.push("C".to_string());
            coords.push([1.39 * a.cos() * b, 1.39 * a.sin() * b, 0.0]);
            rads.push(0.75 * b);
            let a2 = std::f64::consts::TAU * (k as f64 + 0.5) / 6.0;
            syms.push("H".to_string());
            coords.push([2.47 * a2.cos() * b, 2.47 * a2.sin() * b, 0.0]);
            rads.push(0.32 * b);
            bonds.push([2 * k, (2 * k + 2) % 12]);
            bonds.push([2 * k, 2 * k + 1]);
        }
        let mol = Molecule {
            syms,
            coords,
            rads,
            colors: vec![Rgb(142, 142, 142); 12],
            bonds,
            axes: IDENTITY,
        };
        for style in [Style::BallStick, Style::SpaceFilling, Style::Wireframe] {
            let view = View { style, labels: false, ..View::default() };
            let frame = view.render(&mol, 60, 30);
            let w = frame.width as usize;
            let h = frame.height as usize;
            let mut border = Vec::new();
            for col in 0..w {
                border.push(frame.cells[col]);
                border.push(frame.cells[(h - 1) * w + col]);
            }
            for row in 0..h {
                border.push(frame.cells[row * w]);
                border.push(frame.cells[row * w + w - 1]);
            }
            assert!(border.iter().all(|cell| cell.top == BG && cell.bottom == BG), "{:?} 模型把分子画到了画布边缘（自适应缩放没算上球半径）", style);
        }
    }

    #[test]
    fn benzene_from_file() {
        if !Path::new("moles/C6H6.fch").is_file() {
            return;
        }
        let mol = load(Path::new("moles/C6H6.fch")).expect("读取 C6H6 失败");
        assert_eq!(mol.len(), 12);
        assert_eq!(mol.bonds.len(), 12, "苯环应有 12 根键");
        let text = View::default().render(&mol, 80, 30).plain();
        assert!(text.contains('C'));
    }
}
