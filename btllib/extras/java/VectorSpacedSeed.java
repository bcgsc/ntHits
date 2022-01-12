/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package btllib;

public class VectorSpacedSeed extends java.util.AbstractList<VectorUnsigned> implements java.util.RandomAccess {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected VectorSpacedSeed(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(VectorSpacedSeed obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  @SuppressWarnings("deprecation")
  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        btllibJNI.delete_VectorSpacedSeed(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public VectorSpacedSeed(VectorUnsigned[] initialElements) {
    this();
    reserve(initialElements.length);

    for (VectorUnsigned element : initialElements) {
      add(element);
    }
  }

  public VectorSpacedSeed(Iterable<VectorUnsigned> initialElements) {
    this();
    for (VectorUnsigned element : initialElements) {
      add(element);
    }
  }

  public VectorUnsigned get(int index) {
    return doGet(index);
  }

  public VectorUnsigned set(int index, VectorUnsigned e) {
    return doSet(index, e);
  }

  public boolean add(VectorUnsigned e) {
    modCount++;
    doAdd(e);
    return true;
  }

  public void add(int index, VectorUnsigned e) {
    modCount++;
    doAdd(index, e);
  }

  public VectorUnsigned remove(int index) {
    modCount++;
    return doRemove(index);
  }

  protected void removeRange(int fromIndex, int toIndex) {
    modCount++;
    doRemoveRange(fromIndex, toIndex);
  }

  public int size() {
    return doSize();
  }

  public VectorSpacedSeed() {
    this(btllibJNI.new_VectorSpacedSeed__SWIG_0(), true);
  }

  public VectorSpacedSeed(VectorSpacedSeed other) {
    this(btllibJNI.new_VectorSpacedSeed__SWIG_1(VectorSpacedSeed.getCPtr(other), other), true);
  }

  public long capacity() {
    return btllibJNI.VectorSpacedSeed_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    btllibJNI.VectorSpacedSeed_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return btllibJNI.VectorSpacedSeed_isEmpty(swigCPtr, this);
  }

  public void clear() {
    btllibJNI.VectorSpacedSeed_clear(swigCPtr, this);
  }

  public VectorSpacedSeed(int count, VectorUnsigned value) {
    this(btllibJNI.new_VectorSpacedSeed__SWIG_2(count, VectorUnsigned.getCPtr(value), value), true);
  }

  private int doSize() {
    return btllibJNI.VectorSpacedSeed_doSize(swigCPtr, this);
  }

  private void doAdd(VectorUnsigned x) {
    btllibJNI.VectorSpacedSeed_doAdd__SWIG_0(swigCPtr, this, VectorUnsigned.getCPtr(x), x);
  }

  private void doAdd(int index, VectorUnsigned x) {
    btllibJNI.VectorSpacedSeed_doAdd__SWIG_1(swigCPtr, this, index, VectorUnsigned.getCPtr(x), x);
  }

  private VectorUnsigned doRemove(int index) {
    return new VectorUnsigned(btllibJNI.VectorSpacedSeed_doRemove(swigCPtr, this, index), true);
  }

  private VectorUnsigned doGet(int index) {
    return new VectorUnsigned(btllibJNI.VectorSpacedSeed_doGet(swigCPtr, this, index), false);
  }

  private VectorUnsigned doSet(int index, VectorUnsigned val) {
    return new VectorUnsigned(btllibJNI.VectorSpacedSeed_doSet(swigCPtr, this, index, VectorUnsigned.getCPtr(val), val), true);
  }

  private void doRemoveRange(int fromIndex, int toIndex) {
    btllibJNI.VectorSpacedSeed_doRemoveRange(swigCPtr, this, fromIndex, toIndex);
  }

}
