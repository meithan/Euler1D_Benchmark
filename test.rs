
const M: usize = 3;
const N: usize = 2;

fn update_2d_array(arr: &mut Vec<Vec<f64>>) {
  for i in 0..arr.len() {
    for j in 0..arr[i].len() {
      arr[i][j] = 42.0;
    }
  }
}

fn update_passed_by_ref(x: &mut f64) {
  *x = 42.0;
}

fn main() {
  
  let mut array_2d = vec![vec![0f64; M]; N];
  
  println!("Before: {:?}", array_2d);
  update_2d_array(&mut array_2d);
  println!("After: {:?}", array_2d);

  let mut x: f64 = 0.0;
  
  update_passed_by_ref(&mut x);
  println!("x = {}", x);

  let fname = format!("foo{:04}.txt", 123);
  println!("{}", fname);

  for i in 0..5 {
    print!("{} ", i);
  }

  println!();
  for i in 0..=5 {
    print!("{} ", i);
  }

}