use crate::base::Atoms;
use crate::base::Basis;
use crate::base::Shl;
use crate::datas::elems::ELEMS;
use std::collections::HashMap;

// build.rs 会扫描 txts 目录并生成该文件, 因此新增基组文件后无需修改代码
#[cfg(not(feature = "use-rsbase"))]
include!(concat!(env!("OUT_DIR"), "/basis_sets.rs"));

/// 根据基组名和原子数据获取基组对象
///
/// 使用内嵌在 `txts` 目录中的基组数据, 基组名即文件名(不含扩展名),
/// 匹配时忽略大小写和空白字符。只接受基组自身的名字, 不识别别名。
#[cfg(not(feature = "use-rsbase"))]
pub fn from_name(name: &str, atoms: &Atoms) -> Option<Basis> {
    let key = normalize_name(name);
    let text = match BASIS_SETS.iter().find(|(basis, _)| *basis == key) {
        Some((_, text)) => *text,
        None => {
            log::error!("不支持的基组: {}", name);
            return None;
        }
    };
    from_text(text, atoms)
}

/// 根据基组名和原子数据获取基组对象
///
/// 使用 `rsbse`(Basis Set Exchange) 数据库获取基组, 支持其全部基组。
/// 需要先安装 BSE 数据(`bse fetch-data`), 或通过环境变量 `BSE_DATA_DIR` 指定数据目录。
#[cfg(feature = "use-rsbase")]
pub fn from_name(name: &str, atoms: &Atoms) -> Option<Basis> {
    let name = bse_name(name);
    // 必须关闭BSE头部注释, 否则解析器会把它当作未知内容
    let opts = rsbse_ext::Options::new().fmt("gaussian94").header(false);
    let text = match rsbse_ext::get_basis_str(&name, &opts) {
        Ok(text) => text,
        Err(error) => {
            log::error!("获取基组 {} 失败: {}", name, error);
            return None;
        }
    };
    from_text(&text, atoms)
}

/// 当前后端支持的全部基组名（可直接传给 [`from_name`]）
///
/// - 内置后端（未启用 `use-rsbase`）：返回内嵌 `txts` 目录收录的基组名
///   （已归一化：小写、去空白，且按字母序）；
/// - 外部后端（启用 `use-rsbase`）：返回 BSE 数据库中的全部基组名（按字母序）。
#[cfg(not(feature = "use-rsbase"))]
pub fn supported_names() -> Vec<String> {
    BASIS_SETS.iter().map(|(name, _)| (*name).to_string()).collect()
}

/// 当前后端支持的全部基组名（外部 BSE 后端）
#[cfg(feature = "use-rsbase")]
pub fn supported_names() -> Vec<String> {
    match rsbse_ext::get_all_basis_names(None) {
        Ok(names) => names,
        Err(error) => {
            log::error!("获取 BSE 基组列表失败: {}", error);
            Vec::new()
        }
    }
}

/// 把常见写法映射为BSE的基组名
#[cfg(feature = "use-rsbase")]
fn bse_name(name: &str) -> String {
    match normalize_name(name).as_str() {
        "6-31g(d,p)" => "6-31G**".to_string(),
        "6-31g(d)" => "6-31G*".to_string(),
        key => key.to_string(),
    }
}

/// 基组名统一为小写并去掉空白字符
fn normalize_name(name: &str) -> String {
    name.chars().filter(|char| !char.is_whitespace()).collect::<String>().to_lowercase()
}

/// 从基组文件读取基组对象
///
/// 读取 `path` 处的文本后交由 [`from_text`] 解析。
///
/// # Parameters
/// - `path`: 基组文件路径
/// - `atoms`: 分子原子数据，用于筛选并映射基组中的元素
///
/// # Returns
/// 成功时返回与原子数据对应的基组；文件读取或解析失败时返回 `None`
pub fn from_file(path: &str, atoms: &Atoms) -> Option<Basis> {
    let text = std::fs::read_to_string(path).ok()?;
    from_text(&text, atoms)
}

/// 从基组文本解析基组对象
///
/// 文本格式为 Gaussian94 基组格式（元素头 + 壳层头 + 指数/收缩系数数据行），
/// 解析出各元素的壳层后按 `atoms` 中的原子顺序展开，并将元素索引映射为原子索引。
///
/// # Parameters
/// - `text`: 基组文本内容
/// - `atoms`: 分子原子数据，用于筛选基组覆盖的元素并建立映射
///
/// # Returns
/// 成功时返回基组；基组不含分子中某个元素或解析失败时返回 `None`
pub fn from_text(text: &str, atoms: &Atoms) -> Option<Basis> {
    let lines: Vec<String> = text.lines().map(|line| line.to_string()).collect();
    let datas = read_basis(&lines)?;

    // 基组数据覆盖的元素有限, 分子中含有未覆盖的元素时返回None
    let nucs = atoms.nucs();
    for nuc in nucs.iter() {
        let index = nuc.checked_sub(1)?;
        if !datas.iter().any(|data| data.0 == index) {
            log::error!("基组不支持元素 {}", nuc);
            return None;
        }
    }

    // 基组数据按元素记录, 这里按分子中的原子筛选, 并将元素索引映射为原子索引
    let mut shls: Vec<Shl> = Vec::new();
    for (atm, &nuc) in nucs.iter().enumerate() {
        let nuc = nuc - 1;
        let mut keys: Vec<(usize, usize)> = Vec::new(); // (shl, ang)
        let mut atm_shls: Vec<Shl> = Vec::new();
        for &(bnuc, bshl, bang, alp, coe) in datas.iter() {
            if bnuc != nuc {
                continue;
            }
            let key = (bshl, bang);
            if let Some(idx) = keys.iter().position(|&k| k == key) {
                atm_shls[idx].alps.push(alp);
                atm_shls[idx].coes.push(coe);
            } else {
                keys.push(key);
                atm_shls.push(Shl::new(atm, bshl, bang, "car".to_string(), vec![alp], vec![coe]));
            }
        }
        shls.extend(atm_shls);
    }
    Some(Basis::new(shls, None))
}

/// 解析 Gaussian94 格式的基组文本，返回逐条原始数据
///
/// 逐行识别元素头（如 `H 0` 或原子序数）、壳层头（如 `S 3 1.00`、`SP 3 1.00`）
/// 与数据行（指数 + 各角动量对应的收缩系数），最后按（元素、壳层、角动量）排序。
///
/// # Parameters
/// - `lines`: 基组文本按行切分得到的各行字符串
///
/// # Returns
/// 每条数据为 `(元素索引, 壳层序号, 角动量, 指数, 收缩系数)` 的元组列表；
/// 目前总是返回 `Some`
// nuc, shl, ang, alp, coe
pub fn read_basis(lines: &[String]) -> Option<Vec<(usize, usize, usize, f64, f64)>> {
    let mut ang_map: HashMap<char, usize> = HashMap::new();
    for (ang, sym) in ['S', 'P', 'D', 'F', 'G', 'H', 'I'].iter().enumerate() {
        ang_map.insert(*sym, ang);
    }
    let mut datas = Vec::new();
    let mut angs: Vec<usize> = Vec::new();
    let mut shl: usize = 0;
    let mut nuc: usize = 0;
    for line in lines.iter() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('!') || line.starts_with('#') {
            continue;
        }
        // 元素块结束
        if line.starts_with("****") {
            shl = 0;
            angs.clear();
            continue;
        }
        let tokens: Vec<&str> = line.split_whitespace().collect();

        // 元素头, 例如 "H 0" (Gaussian基组文件) 或 "1 0" (Gaussian日志)
        if tokens.len() == 2 {
            if let Ok(flag) = tokens[1].parse::<i64>() {
                if let Ok(atm) = tokens[0].parse::<usize>() {
                    // 日志中以原子序数开头
                    nuc = atm - 1;
                    shl = 0;
                    angs.clear();
                    continue;
                }
                // 单字母元素符号与角动量符号(S/P/F)有歧义, 因此要求第二个数为0
                if let Some(index) = elem_index(tokens[0]) {
                    if flag == 0 || tokens[0].len() > 1 {
                        nuc = index;
                        shl = 0;
                        angs.clear();
                        continue;
                    }
                }
            }
        }

        // 壳层头, 例如 "S 3 1.00" 或 "SP 3 1.00", 日志中末尾还带一列零
        if tokens.len() >= 2 && tokens[1].parse::<usize>().is_ok() {
            let shell = tokens[0].to_ascii_uppercase();
            if !shell.is_empty() && shell.chars().all(|sym| ang_map.contains_key(&sym)) {
                angs = shell.chars().map(|sym| ang_map[&sym]).collect();
                shl += 1;
                continue;
            }
        }

        // 数据行: 指数 + 各角动量对应的收缩系数
        if tokens.len() >= 2 {
            let nums: Option<Vec<f64>> = tokens.iter().map(|each| parse_float(each)).collect();
            if let Some(nums) = nums {
                if !angs.is_empty() && shl > 0 {
                    let alp = nums[0];
                    for (idx, ang) in angs.iter().enumerate() {
                        if let Some(&coe) = nums.get(idx + 1) {
                            datas.push((nuc, shl - 1, *ang, alp, coe));
                        }
                    }
                }
                continue;
            }
        }

        // 无法识别的内容, 认为基组数据已结束
        break;
    }
    datas.sort_by_key(|data| (data.0, data.1, data.2));
    Some(datas)
}

/// 根据元素符号返回从0开始的元素索引
fn elem_index(sym: &str) -> Option<usize> {
    ELEMS.iter().find(|elm| elm.idx > 0 && elm.sym.eq_ignore_ascii_case(sym)).map(|elm| elm.idx as usize - 1)
}

/// 解析基组中的浮点数, 兼容Fortran的D指数
fn parse_float(text: &str) -> Option<f64> {
    text.replace('D', "e").replace('d', "e").parse::<f64>().ok()
}

#[cfg(test)]
mod tests {
    use super::{from_name, read_basis};
    use crate::base::Atoms;
    use std::collections::BTreeSet;
    use std::fs;

    #[test]
    fn test_read_gaussian_basis() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/src/rsbse/txts/6-31g(d,p).txt");
        let content = fs::read_to_string(path).expect("读取基组文件失败");
        let lines: Vec<String> = content.lines().map(|line| line.to_string()).collect();
        let datas = read_basis(&lines).expect("解析基组失败");
        assert!(!datas.is_empty());

        // H 是最小索引0, 共 3+1 个s和1个p
        let h: Vec<_> = datas.iter().filter(|data| data.0 == 0).collect();
        assert_eq!(h.len(), 5);
        assert!(h.iter().any(|data| data.2 == 1));
        assert!((h[0].3 - 18.73113696).abs() < 1e-10);
        assert!((h[0].4 - 0.03349460434).abs() < 1e-10);

        // He 的索引为1
        assert!(datas.iter().any(|data| data.0 == 1));

        // Li 的索引为2, 含SP壳层
        let li: Vec<_> = datas.iter().filter(|data| data.0 == 2).collect();
        assert!(li.iter().any(|data| data.2 == 0));
        assert!(li.iter().any(|data| data.2 == 1));

        // 文件覆盖 H-Kr 共36个元素, 每个指数在每个角动量上各有一个系数
        let nucs: BTreeSet<usize> = datas.iter().map(|data| data.0).collect();
        assert_eq!(nucs.len(), 36);
        assert_eq!(nucs.iter().max(), Some(&35));
        assert_eq!(datas.len(), 1120);
    }

    #[cfg(not(feature = "use-rsbase"))]
    #[test]
    fn test_from_name() {
        // 水分子: O, H, H
        let atms = vec![8, 1, 1];
        let xyzs = vec![[0.0, 0.0, 0.0], [0.0, 0.757, 0.586], [0.0, -0.757, 0.586]];
        let crgs = vec![8.0, 1.0, 1.0];
        let atoms = Atoms::new(&atms, &xyzs, &crgs);

        let basis = from_name("6-31G(d,p)", &atoms).expect("获取基组失败");
        // O 有6个壳层(S, SP->S/P, SP->S/P, D), 每个H有3个壳层
        assert_eq!(basis.shls.len(), 12);
        // 水分子在笛卡尔形式下共 15(O) + 5*2(H) = 25 个基函数
        assert_eq!(basis.shls.nato("car"), 25);
        assert_eq!(basis.shls.nato("sph"), 24);
        for atm in 0..atoms.len() {
            assert!(basis.shls.shls.iter().any(|shl| shl.atm == atm));
        }
        // 所有壳层使用笛卡尔形式
        assert!(basis.shls.shls.iter().all(|shl| shl.form == "car"));

        // 只接受基组文件自身的名字, 不识别别名
        assert!(from_name("6-31g**", &atoms).is_none());
        // txts 下的文件会被自动收录, 无需修改代码
        let def2 = from_name("def2-SVP", &atoms).expect("获取def2-SVP基组失败");
        for atm in 0..atoms.len() {
            assert!(def2.shls.shls.iter().any(|shl| shl.atm == atm));
        }
        // 不存在的基组返回None
        assert!(from_name("cc-pvdz", &atoms).is_none());
        // 含未覆盖元素(如 Xe)时返回None
        let atms = vec![54];
        let xyzs = vec![[0.0, 0.0, 0.0]];
        let crgs = vec![54.0];
        let atoms = Atoms::new(&atms, &xyzs, &crgs);
        assert!(from_name("6-31G(d,p)", &atoms).is_none());
    }

    #[cfg(not(feature = "use-rsbase"))]
    #[test]
    fn test_supported_names() {
        let names = super::supported_names();
        assert!(!names.is_empty());
        // 列表与 txts 目录一致（归一化后的名字）
        for expect in ["6-31g(d,p)", "def2-svp", "sto-3g"] {
            assert!(names.iter().any(|name| name.as_str() == expect), "缺少基组 {expect}");
        }

        // 每个名字都能被 from_name 解析（水分子覆盖 O、H）
        let atms = vec![8, 1, 1];
        let xyzs = vec![[0.0, 0.0, 0.0], [0.0, 0.757, 0.586], [0.0, -0.757, 0.586]];
        let crgs = vec![8.0, 1.0, 1.0];
        let atoms = Atoms::new(&atms, &xyzs, &crgs);
        for name in &names {
            assert!(from_name(name, &atoms).is_some(), "基组 {name} 无法解析");
        }
    }

    #[cfg(feature = "use-rsbase")]
    #[test]
    fn test_from_name_bse() {
        // 需要先安装BSE数据: `bse fetch-data`, 或设置环境变量 BSE_DATA_DIR
        let data_dir = rsbse_ext::get_data_dir();
        if !data_dir.join("METADATA.json").is_file() {
            eprintln!("跳过测试: 未安装BSE数据 ({})", data_dir.display());
            return;
        }

        // 水分子: O, H, H
        let atms = vec![8, 1, 1];
        let xyzs = vec![[0.0, 0.0, 0.0], [0.0, 0.757, 0.586], [0.0, -0.757, 0.586]];
        let crgs = vec![8.0, 1.0, 1.0];
        let atoms = Atoms::new(&atms, &xyzs, &crgs);

        let basis = from_name("6-31G(d,p)", &atoms).expect("获取基组失败");
        // 与内嵌数据的结果一致
        assert_eq!(basis.shls.len(), 12);
        assert_eq!(basis.shls.nato("car"), 25);
        for atm in 0..atoms.len() {
            assert!(basis.shls.shls.iter().any(|shl| shl.atm == atm));
        }
        // BSE 数据库中使用 * 表示极化函数
        assert!(from_name("6-31G**", &atoms).is_some());
        assert!(from_name("def2-SVP", &atoms).is_some());
        assert!(from_name("sto-3g", &atoms).is_some());
        // 不存在的基组返回None
        assert!(from_name("not-a-basis", &atoms).is_none());
    }

    #[cfg(feature = "use-rsbase")]
    #[test]
    fn test_supported_names_bse() {
        let data_dir = rsbse_ext::get_data_dir();
        if !data_dir.join("METADATA.json").is_file() {
            eprintln!("跳过测试: 未安装BSE数据 ({})", data_dir.display());
            return;
        }

        let names = super::supported_names();
        assert!(!names.is_empty());
        // BSE 返回显示名（大小写不敏感比较）
        assert!(names.iter().any(|name| name.eq_ignore_ascii_case("def2-SVP")));
        assert!(names.iter().any(|name| name.eq_ignore_ascii_case("6-31G*")));
    }
}
