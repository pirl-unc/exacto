use crate::structs::trie::Trie;


#[test]
fn test_trie_1() {
    let mut trie: Trie = Trie::new();
    trie.insert("apple");
    trie.insert("appetizer");
    trie.insert("application");
    assert!(trie.exists("app"));
    assert!(trie.exists("appl"));
    assert!(trie.exists("apple"));
    assert!(trie.search("app").len() == 3);
}
