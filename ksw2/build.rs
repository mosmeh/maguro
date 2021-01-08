fn main() {
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
}
