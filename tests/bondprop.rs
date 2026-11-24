use rswfn::base::Mole;
use rswfn::bondprop;

mod utils;
use utils::allclose;

#[test]
/// 测试含金属环的log文件计算的mayer键级是否正确
fn test1() {
    let path = "moles/dui-Nb.log_scan_14.log";
    let mole = Mole::from_file(&path);
    let caler = bondprop::order::Calculator::new(&mole);
    let omat = caler.mayer();
    let natm = mole.atoms.len();
    let mut res = vec![];
    for i in 0..natm - 1 {
        for j in i + 1..natm {
            let val = omat[(i, j)];
            if val > 0.05 {
                res.push(val);
            }
        }
    }
    #[rustfmt::skip]
    let ans=vec![
        1.35469760,1.45377985,0.94309473,0.15647486,1.35460700,1.05424435,0.12685757,1.45404800,0.85518334,1.44112431,
        0.94309644,0.15647474,0.08199164,0.85506430,1.44158642,0.08217649,1.35462738,1.45394850,0.94309817,0.15647280,
        1.35467197,0.12686021,1.45393458,0.85514935,1.44130197,0.94309388,0.15646002,0.08204872,0.85512087,1.44132495,
        0.08208664,0.08512342,
    ];
    assert!(allclose(res, ans, 1e-4));
}
