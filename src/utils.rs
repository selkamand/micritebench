pub fn prefix(path: &str) -> String {
    let mypath = std::path::Path::new(&path);
    let filename = match mypath.file_name() {
        Some(val) => val,
        None => panic!("Failed to get file prefix"),
    };

    filename
        .to_string_lossy()
        .split('.')
        .next()
        .unwrap_or("")
        .to_string()
}
