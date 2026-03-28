use std::collections::HashMap;

fn main() {
    let mut hm = HashMap::new();
    hm.insert(b"read1".to_vec(), 0);
    println!("{:?}", hm);
}
