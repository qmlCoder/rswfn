use std::f64::consts::PI;
use std::fmt;
use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::MulAssign;
use std::ops::Sub;

/// 三维向量，用于原子坐标、方向向量等运算
#[derive(Debug, Clone)]
pub struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3 {
    /// 创建一个新的三维向量
    ///
    /// # Parameters
    /// - `x`: x 分量
    /// - `y`: y 分量
    /// - `z`: z 分量
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// 从长度为3的数组创建三维向量
    ///
    /// 会检查所有分量是否为非 NaN 值
    ///
    /// # Parameters
    /// - `vec`: [x, y, z] 数组
    pub fn from_arr(vec: [f64; 3]) -> Self {
        let x = vec[0];
        let y = vec[1];
        let z = vec[2];
        assert!(!x.is_nan() && !y.is_nan() && !z.is_nan());
        Vec3::new(x, y, z)
    }

    /// 从 `Vec<f64>` 创建三维向量
    ///
    /// 会检查所有分量是否为非 NaN 值
    ///
    /// # Parameters
    /// - `vec`: 至少包含3个元素的向量
    pub fn from_vec(vec: &Vec<f64>) -> Self {
        let x = vec[0];
        let y = vec[1];
        let z = vec[2];
        assert!(!x.is_nan() && !y.is_nan() && !z.is_nan());
        Vec3::new(x, y, z)
    }

    /// 将向量转换为长度为3的数组
    pub fn into_arr(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// 将向量转换为 `Vec<f64>`
    pub fn into_vec(&self) -> Vec<f64> {
        vec![self.x, self.y, self.z]
    }

    /// 返回零向量 (0, 0, 0)
    pub fn zero() -> Self {
        Vec3::new(0.0, 0.0, 0.0)
    }

    /// 生成一个随机单位向量
    ///
    /// 每个分量在 [0, 1) 范围内随机生成，然后归一化
    pub fn rand() -> Self {
        let x: f64 = rand::random();
        let y: f64 = rand::random();
        let z: f64 = rand::random();
        let l = (x * x + y * y + z * z).sqrt();
        Vec3::new(x / l, y / l, z / l)
    }

    /// 返回归一化后的单位向量
    pub fn normalize(&self) -> Vec3 {
        let norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        Vec3::new(self.x / norm, self.y / norm, self.z / norm)
    }

    /// 返回向量的长度（模）
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// 计算与另一个向量的叉积
    ///
    /// # Parameters
    /// - `vec3`: 另一个三维向量
    pub fn cross(&self, vec3: &Vec3) -> Vec3 {
        let [x1, y1, z1] = self.into_arr();
        let [x2, y2, z2] = vec3.into_arr();
        Vec3::new(y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2)
    }

    /// 计算与另一个向量的点积
    ///
    /// # Parameters
    /// - `vec3`: 另一个三维向量
    pub fn dot(&self, vec3: &Vec3) -> f64 {
        self.x * vec3.x + self.y * vec3.y + self.z * vec3.z
    }

    /// 两向量之间的角度 0-PI
    pub fn angle(&self, vec3: &Vec3) -> f64 {
        let cos = self.dot(vec3) / (self.length() * vec3.length());
        let cos = cos.max(-1.0).min(1.0);
        cos.acos()
    }

    /// 计算向量在另一个向量方向上的投影
    ///
    /// # Parameters
    /// - `vec3`: 目标投影方向向量
    pub fn project_onto(&self, vec3: &Vec3) -> Vec3 {
        vec3 * self.dot(vec3) / vec3.length()
    }

    /// 旋转一个点，self代表point
    pub fn rotate_arrow(&self, axis: &Vec3, angle: f64) -> Self {
        let axis = axis.normalize();
        let [px, py, pz] = self.into_arr();
        let [dx, dy, dz] = axis.into_arr();
        let cos = angle.cos();
        let sin = angle.sin();
        let x = px * (dx.powi(2) * (1.0 - cos) + cos)
            + py * (dx * dy * (1.0 - cos) - dz * sin)
            + pz * (dx * dz * (1.0 - cos) + dy * sin);
        let y = px * (dx * dy * (1.0 - cos) + dz * sin)
            + py * (dy.powi(2) * (1.0 - cos) + cos)
            + pz * (dy * dz * (1.0 - cos) - dx * sin);
        let z = px * (dx * dz * (1.0 - cos) - dy * sin)
            + py * (dy * dz * (1.0 - cos) + dx * sin)
            + pz * (dz.powi(2) * (1.0 - cos) + cos);
        Vec3::new(x, y, z)
    }

    /// 旋转箭头，self代表箭头的朝向
    pub fn rotate_point(&self, cent: &Vec3, axis: &Vec3, angle: f64) -> Self {
        let point = self - cent;
        let point = point.rotate_arrow(axis, angle);
        point + cent
    }

    /// 计算点到平面的距离
    ///
    /// 平面由方程 ax + by + cz + d = 0 定义
    ///
    /// # Parameters
    /// - `plan`: 平面参数 [a, b, c, d]
    pub fn dist_to_plan(&self, plan: &[f64; 4]) -> f64 {
        let [a, b, c, d] = plan;
        (a * self.x + b * self.y + c * self.z + d) / (a * a + b * b + c * c).sqrt()
    }
}

impl Sub<&Vec3> for &Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Sub<Vec3> for Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Sub<&Vec3> for Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Sub<Vec3> for &Vec3 {
    type Output = Vec3;
    fn sub(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x - vec3.x, self.y - vec3.y, self.z - vec3.z)
    }
}

impl Add<&Vec3> for &Vec3 {
    type Output = Vec3;
    fn add(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Add<Vec3> for Vec3 {
    type Output = Vec3;
    fn add(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Add<&Vec3> for Vec3 {
    type Output = Vec3;
    fn add(self, vec3: &Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Add<Vec3> for &Vec3 {
    type Output = Vec3;
    fn add(self, vec3: Vec3) -> Vec3 {
        Vec3::new(self.x + vec3.x, self.y + vec3.y, self.z + vec3.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, scalar: f64) -> Vec3 {
        Vec3::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Mul<f64> for &Vec3 {
    type Output = Vec3;
    fn mul(self, scalar: f64) -> Vec3 {
        Vec3::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, scalar: f64) -> Vec3 {
        Vec3::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }
}

impl MulAssign<f64> for Vec3 {
    fn mul_assign(&mut self, scalar: f64) {
        self.x *= scalar;
        self.y *= scalar;
        self.z *= scalar;
    }
}

impl fmt::Display for Vec3 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "({:>10.4}, {:>10.4}, {:>10.4}):{:>10.4}",
            self.x,
            self.y,
            self.z,
            self.length()
        )?;
        Ok(())
    }
}
