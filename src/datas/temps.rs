/// Gaussian 输入文件（.gjf）模板
pub static GJF: &str = r#"%chk={chk}.chk
# {job}

Generate by Pywfn

{ele.0} {ele.1}
{geome}


"#;

/// XYZ 坐标文件模板
pub static XYZ: &str = r#"{natm}
{title}
{geome}
"#;

/// 补充材料（SI）文本模板（坐标、能量与振动频率等）
pub static SI: &str = r#"{name}
{engs}

Cartesian coordinate
{geome}

Variable frequencies
{freqs}
"#;
