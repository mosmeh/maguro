use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=wrapper.h");

    cc::Build::new()
        .warnings(false)
        .file("ksw2/kalloc.c")
        .file("ksw2/ksw2_gg.c")
        .file("ksw2/ksw2_gg2.c")
        .file("ksw2/ksw2_gg2_sse.c")
        .file("ksw2/ksw2_extz.c")
        .file("ksw2/ksw2_extz2_sse.c")
        .file("ksw2/ksw2_extd.c")
        .file("ksw2/ksw2_extd2_sse.c")
        .compile("ksw2");

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
