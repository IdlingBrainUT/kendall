use avlsort::tree::*;
use avlsort::traits::*;

use std::cmp::Ordering;
use std::collections::BinaryHeap;

pub fn tau_def(c: i128, d: i128, ex: i128, ey: i128) -> f64 {
    let (cf, df, exf, eyf) = (c as f64, d as f64, ex as f64, ey as f64);
    (cf - df) / ((cf + df + exf).sqrt() * (cf + df + eyf).sqrt())
}

pub fn tau_scale_def(n: f64, c: i128, d: i128, ex: i128, ey: i128, ex1: f64, ey1: f64, ex2: f64, ey2: f64) -> f64 {
    if n <= 2.0 {
        panic!();
    }
    let n_n1 = n * (n - 1.0);
    let n_n1_n2 = n_n1 * (n - 2.0);
    let n_n1_2n_5 = n_n1 * (2.0 * n + 5.0);
    let var = (n_n1_2n_5 - ex2 - ey2) / 18.0 
                + (ex1 * ey1) / (9.0 * n_n1_n2) 
                + (ex * ey) as f64 / (2.0 * n_n1);
    (c - d) as f64 / var.sqrt()
}

pub fn tau<T: TreeElem>(a: &[T], b: &[T]) -> Result<f64, ()> {
    let (c, d, ex, ey) = tau_core(a, b)?;
    Ok(tau_def(c, d, ex, ey))
}

pub fn tau_core<T: TreeElem>(a: &[T], b: &[T]) -> Result<(i128, i128, i128, i128), ()> {
    if a.len() != b.len() {
        return Err(());
    }
    let (mut c, mut d, mut ex, mut ey) = (0, 0, 0, 0);
    for (i, (&ai, &bi)) in a.iter().zip(b.iter()).enumerate() {
        for (&aj, &bj) in a[i+1..].iter().zip(b[i+1..].iter()) {
            if ai < aj {
                if bi < bj {
                    c += 1;
                } else if bi > bj {
                    d += 1;
                } else {
                    ex += 1;
                }
            } else if ai > aj {
                if bi < bj {
                    d += 1;
                } else if bi > bj {
                    c += 1;
                } else {
                    ex += 1;
                }
            } else {
                if bi != bj {
                    ey += 1;
                }
            }
        }
    }

    Ok((c, d, ex, ey))
}

pub fn ndtau<T: TreeElem>(a: &[T], b: &[T]) -> Result<f64, ()> {
    let (n, c) = ndtau_core(a, b)?;
    let (nf, cf) = (n as f64, c as f64);
    Ok(4.0 * cf / (nf * (nf - 1.0)) - 1.0)
}

pub fn ndtau_core<T: TreeElem>(a: &[T], b: &[T]) -> Result<(i128, i128), ()> {
    let n = a.len() as i128;
    if n != b.len() as i128 {
        return Err(());
    }

    let mut ab = a.iter().zip(b.iter()).map(|(&ai, &bi)| (ai, bi)).collect::<Vec<(T, T)>>();
    ab.sort_by(|s, t| ord_tup(s, t));
    let mut avl = AvlTree::new();
    let mut c = 0i128;
    for &(_, bi) in ab.iter() {
        c += avl.push(bi).0 as i128;
    }
    Ok((n, c))
}

pub fn sdtau<T: TreeElem>(a: &[T], b: &[T]) -> Result<f64, ()> {
    let (c, d, ex, ey) = sdtau_core(a, b)?;
    Ok(tau_def(c, d, ex, ey))
}

pub fn sdtau_core<T: TreeElem>(a: &[T], b: &[T]) -> Result<(i128, i128, i128, i128), ()> {
    let n = a.len();
    if n != b.len() {
        return Err(());
    }

    let mut ab = a.iter().zip(b.iter()).map(|(&ai, &bi)| (ai, bi)).collect::<Vec<(T, T)>>();
    ab.sort_by(|s, t| ord_tup(s, t));
    let mut avl = AvlTree::new();
    let (mut c, mut d, mut ex, mut ey) = (0i128, 0i128, 0i128, 0i128);
    let (mut d_count, mut e_count) = (0i128, 1i128);
    let mut prev = ab[0];
    for (i, &abi) in ab.iter().enumerate() {
        if i > 0 {
            if abi.0 != prev.0 {
                d_count = 0;
                e_count = 1;
            } else {
                if abi.1 == prev.1 {
                    e_count += 1;
                } else {
                    d_count += e_count;
                    e_count = 1;
                }
            }
        }
        prev = abi;
        let (nb, ne) = avl.push(abi.1);
        let a_count = nb as i128 - d_count;
        let b_count = (ne + 1) as i128 - e_count;
        let c_count = (i + 1) as i128 - (a_count + b_count + d_count + e_count);
        ey += d_count;
        ex += b_count;
        c += a_count;
        d += c_count;
    }

    Ok((c, d, ex, ey))
}

pub fn nztau<T: TreeElem>(a: &[T], b: &[T]) -> Result<f64, ()> {
    let (c, d, ex, ey) = nztau_core(a, b)?;
    Ok(tau_def(c, d, ex, ey))
}

pub fn nztau_core<T: TreeElem>(a: &[T], b: &[T]) -> Result<(i128, i128, i128, i128), ()> {
    let n = a.len();
    if n != b.len() {
        return Err(());
    }
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut xz = AvlTree::new();
    let mut zy = AvlTree::new();
    let mut z = 0;
    let zero = T::zero();

    for (&ai, &bi) in a.iter().zip(b.iter()) {
        if ai > zero {
            if bi > zero {
                x.push(ai);
                y.push(bi);
            } else if bi == zero {
                xz.push(ai);
            } else {
                return Err(());
            }
        } else if ai == zero {
            if bi > zero {
                zy.push(bi);
            } else if bi == zero {
                z += 1;
            } else {
                return Err(());
            }
        } else {
            return Err(());
        }
    }

    let (mut c, mut d, mut ex, mut ey) = sdtau_core(&x, &y)?;

    x.sort_by(|&a, b| a.partial_cmp(b).unwrap());
    let x_len = x.len();
    let xz_len = xz.len();
    let mut xz_comb = if xz_len == 0 { 0 } else { xz_len * (xz_len - 1) / 2 };
    let mut count_less = 0;
    let mut count_same = 0;
    loop {
        match xz.pop_min_all() {
            Some((x0, x0_dup)) => {
                let x0_num = x0_dup + 1;
                xz_comb -= x0_num * (x0_num - 1) / 2;
                loop {
                    if count_less >= x_len || x0 <= x[count_less] { break; }
                    count_same = 0;
                    count_less += 1;
                }
                loop {
                    if count_less + count_same >= x_len || x0 < x[count_less + count_same] { break; }
                    count_same += 1;
                }
                d += (count_less * x0_num) as i128;
                c += ((x_len - count_less - count_same) * x0_num) as i128;
                ey += (count_same * x0_num) as i128;
            }
            None => break,
        }
    }
    ex += xz_comb as i128;

    y.sort_by(|&a, b| a.partial_cmp(b).unwrap());
    let y_len = y.len();
    let zy_len = zy.len();
    let mut zy_comb = if zy_len == 0 { 0 } else { zy_len * (zy_len - 1) / 2 };
    let mut count_less = 0;
    let mut count_same = 0;
    loop {
        match zy.pop_min_all() {
            Some((y0, y0_dup)) => {
                let y0_num = y0_dup + 1;
                zy_comb -= y0_num * (y0_num - 1) / 2;
                loop {
                    if count_less >= y_len || y0 <= y[count_less] { break; }
                    count_same = 0;
                    count_less += 1;
                }
                loop {
                    if count_less + count_same >= y_len || y0 < y[count_less + count_same] { break; }
                    count_same += 1;
                }
                d += (count_less * y0_num) as i128;
                c += ((y_len - count_less - count_same) * y0_num) as i128;
                ex += (count_same * y0_num) as i128;
            }
            None => break,
        }
    }
    ey += zy_comb as i128;

    c += (z * x_len) as i128;

    d += (xz_len * zy_len) as i128;

    ex += (z * xz_len) as i128;

    ey += (z * zy_len) as i128;
    
    Ok((c, d, ex, ey))
}

pub fn ndnztau<T: TreeElem>(a: &[T], b: &[T]) -> Result<f64, ()> {
    let (c, d, ex, ey) = nzndtau_core(a, b)?;
    Ok(tau_def(c, d, ex, ey))
}

pub fn ndnztau_scale<T: TreeElem>(a: &[T], b: &[T]) -> Result<f64, ()> {
    let (c, d, ex, ey) = nzndtau_core(a, b)?;
    let n = a.len() as f64;
    let (ex1, ex2) = if ex == 0 {
        (0.0, 0.0)
    } else {
        let exf = ex as f64;
        let exf_sqrt = exf.sqrt();
        let x = exf_sqrt.ceil();
        let x_x1 = x * (x - 1.0);
        (x_x1 * (x - 2.0), x_x1 * (2.0 * x + 5.0))
    };
    let (ey1, ey2) = if ey == 0 {
        (0.0, 0.0)
    } else {
        let eyf = ey as f64;
        let eyf_sqrt = eyf.sqrt();
        let y = eyf_sqrt.ceil();
        let y_y1 = y * (y - 1.0);
        (y_y1 * (y - 2.0), y_y1 * (2.0 * y + 5.0))
    };
    Ok(tau_scale_def(n, c, d, ex, ey, ex1, ey1, ex2, ey2))
}

pub fn nzndtau_core<T: TreeElem>(a: &[T], b: &[T]) -> Result<(i128, i128, i128, i128), ()> {
    let n = a.len();
    if n != b.len() {
        return Err(());
    }
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut xz = BinaryHeap::new();
    let mut zy = BinaryHeap::new();
    let mut z = 0;
    let zero = T::zero();

    for (&ai, &bi) in a.iter().zip(b.iter()) {
        if ai > zero {
            if bi > zero {
                x.push(OrdEqElem::new(ai));
                y.push(OrdEqElem::new(bi));
            } else if bi == zero {
                xz.push(OrdEqElem::new(ai));
            } else {
                return Err(());
            }
        } else if ai == zero {
            if bi > zero {
                zy.push(OrdEqElem::new(bi));
            } else if bi == zero {
                z += 1;
            } else {
                return Err(());
            }
        } else {
            return Err(());
        }
    }

    let (n, mut c) = ndtau_core(&x, &y)?;
    let mut d = if n <= 1 {
        panic!()
    } else if n % 2 == 0 {
        n / 2 * (n - 1) - c
    } else {
        (n - 1) / 2 * n - c
    };
    let (mut ex, mut ey) = (0i128, 0i128);

    x.sort();
    let x_len = x.len();
    let xz_len = xz.len();
    let xz_comb = if xz_len == 0 { 0 } else { xz_len * (xz_len - 1) / 2 };
    let mut count_less = x_len;
    loop {
        match xz.pop() {
            Some(x0) => {
                loop {
                    if count_less == 0 || x0 > x[count_less - 1] { break; }
                    count_less -= 1;
                }
                d += count_less as i128;
                c += (x_len - count_less) as i128;
            }
            None => break,
        }
    }
    ex += xz_comb as i128;

    y.sort();
    let y_len = y.len();
    let zy_len = zy.len();
    let zy_comb = if zy_len == 0 { 0 } else { zy_len * (zy_len - 1) / 2 };
    let mut count_less = y_len;
    loop {
        match zy.pop() {
            Some(y0) => {
                loop {
                    if count_less == 0 || y0 > y[count_less - 1] { break; }
                    count_less -= 1;
                }
                d += count_less as i128;
                c += (y_len - count_less) as i128;
            }
            None => break,
        }
    }
    ey += zy_comb as i128;

    c += (z * x_len) as i128;

    d += (xz_len * zy_len) as i128;

    ex += (z * xz_len) as i128;

    ey += (z * zy_len) as i128;
    
    Ok((c, d, ex, ey))
}

pub fn ord_tup<T: TreeElem>(s: &(T, T), t: &(T, T)) -> Ordering {
    if s.0 < t.0 {
        Ordering::Less
    } else if s.0 > t.0 {
        Ordering::Greater
    } else {
        if s.1 < t.1 {
            Ordering::Less
        } else if s.1 > t.1 {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}