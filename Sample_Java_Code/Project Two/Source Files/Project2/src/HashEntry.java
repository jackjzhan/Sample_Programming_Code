
public class HashEntry {
	private int key;
    private HashEntry next;

    HashEntry(int key) {
          this.key = key;
          this.next = null;
    }

    public int getKey() {
          return key;
    }

    public HashEntry getNext() {
          return next;
    }

    public void setNext(HashEntry next) {
          this.next = next;
    }
}
