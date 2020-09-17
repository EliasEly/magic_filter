require 'BOAST'
require 'narray_ffi'
include BOAST

# Redefine a simple c modulo
class Modulo
  def to_s_c
    op1, op2 = op_to_var
    return "((#{op1}) % (#{op2}))"
  end
end

CFLAGS = "-O3 -march=native"

FILTER_SIZE = 16

FILTER = [
   8.4334247333529341094733325815816e-7,
  -0.1290557201342060969516786758559028e-4,
   0.8762984476210559564689161894116397e-4,
  -0.30158038132690463167163703826169879e-3,
   0.174723713672993903449447812749852942e-2,
  -0.942047030201080385922711540948195075e-2,
   0.2373821463724942397566389712597274535e-1,
   0.612625895831207982195380597e-1,
   0.9940415697834003993178616713,
  -0.604895289196983516002834636e-1,
  -0.2103025160930381434955489412839065067e-1,
   0.1337263414854794752733423467013220997e-1,
  -0.344128144493493857280881509686821861e-2,
   0.49443227688689919192282259476750972e-3,
  -0.5185986881173432922848639136911487e-4,
   2.72734492911979659657715313017228e-6 ]

FILTER_REVERSE = [
   2.72734492911979659657715313017228e-6,
  -0.5185986881173432922848639136911487e-4,
   0.49443227688689919192282259476750972e-3,
  -0.344128144493493857280881509686821861e-2,
   0.1337263414854794752733423467013220997e-1,
  -0.2103025160930381434955489412839065067e-1,
  -0.604895289196983516002834636e-1,
   0.9940415697834003993178616713,
   0.612625895831207982195380597e-1,
   0.2373821463724942397566389712597274535e-1,
  -0.942047030201080385922711540948195075e-2,
   0.174723713672993903449447812749852942e-2,
  -0.30158038132690463167163703826169879e-3,
   0.8762984476210559564689161894116397e-4,
  -0.1290557201342060969516786758559028e-4,
   8.4334247333529341094733325815816e-7 ]


src_filters = <<EOF
const double filter[] __attribute__ ((aligned (16))) = { #{FILTER.join(", ")} };
const double filter_reverse[] __attribute__ ((aligned (16))) = { #{FILTER_REVERSE.join(", ")} };
EOF

src_magicfilter1d_ref = src_filters.dup
src_magicfilter1d_ref << <<EOF
void magicfilter1d_naive_(int *n, int *ndat, double const *source, double *dest) {
  double tmp;
  unsigned int i,j,k;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      tmp=0;
      for(k=0;k<16;k++){
        tmp+=source[(j-8+k+(*n))%(*n)]*filter[k];
      }
      dest[j*(*ndat)]=tmp;
    }
    dest += 1;
    source += (*n);
  }
}
EOF

src_magicfilter1d_t_ref = src_filters.dup
src_magicfilter1d_t_ref << <<EOF
void magicfilter1d_t_naive_(unsigned int *n, unsigned int *ndat, double const *source, double *dest) {
  double tmp;
  unsigned int i,j,k;
  for(i=0;i<(*ndat);i++){
    for(j=0;j<(*n);j++) {
      tmp=0;
      for(k=0;k<16;k++){
        tmp+=source[(j-7+k+(*n))%(*n)]*filter_reverse[k];
      }
      dest[j*(*ndat)]=tmp;
    }

    dest += 1;
    source += (*n);
  }
}
EOF

set_lang(C)
set_array_start(0)
set_default_real_size(8)

n    = Int :n,    reference: true, dir: :in
ndat = Int :ndat, reference: true, dir: :in
source = Real :source, dim: [Dim(n), Dim(ndat)], dir: :in
dest   = Real :dest,   dim: [Dim(ndat), Dim(n)], dir: :out
magicfilter1d_naive = Procedure(:magicfilter1d_naive_, [n, ndat, source, dest])
magicfilter1d_ref = CKernel::new {
  get_output.print src_magicfilter1d_ref
}
magicfilter1d_ref.procedure = magicfilter1d_naive
magicfilter1d_ref.build(CFLAGS: CFLAGS)

magicfilter1d_t_naive = Procedure(:magicfilter1d_t_naive_, [n, ndat, source, dest])
magicfilter1d_t_ref = CKernel::new {
  get_output.print src_magicfilter1d_t_ref
}
magicfilter1d_t_ref.procedure = magicfilter1d_t_naive
magicfilter1d_t_ref.build(CFLAGS: CFLAGS)

X = 31
Y = 32
Z = 33
TOTAL = (2*X)*(2*Y)*(2*Z)

source = NArray.float(TOTAL).random!
tmp = NArray.float(TOTAL)
dest_ref = NArray.float(TOTAL)
dest_t_ref = NArray.float(TOTAL)
dest = NArray.float(TOTAL)

def performance(res, x, y, z)
  TOTAL * FILTER_SIZE*2 / (res[:duration] * 1e9)
end

def apply_kernel(k, x, y, z, src, tmp, dest)
  res = []
  res.push k.run(2*x, 4*y*z, src, dest)
  res.push k.run(2*y, 4*z*x, dest, tmp)
  res.push k.run(2*z, 4*x*y, tmp, dest)
  res.collect { |r| performance(r, x, y, z) }.reduce(:+) / 3.0
end

NUM_REPEAT = 10

def bench_kernel(k, x, y, z, src, tmp, dest)
  res = NUM_REPEAT.times.collect {
    apply_kernel(k, x, y, z, src, tmp, dest)
  }
  res.max
end

puts "magicfilter1d_naive_ #{X} #{Y} #{Z}"
apply_kernel(magicfilter1d_ref, X, Y, Z, source, tmp, dest_ref)
puts "#{bench_kernel(magicfilter1d_ref, X, Y, Z, source, tmp, dest)} GFlop/s"
puts "magicfilter1d_t_naive_ #{X} #{Y} #{Z}"
apply_kernel(magicfilter1d_t_ref, X, Y, Z, source, tmp, dest_t_ref)
puts "#{bench_kernel(magicfilter1d_t_ref, X, Y, Z, source, tmp, dest)} GFlop/s"

def print_inner_loop(i, j, k, n, tmp, offset, source, dest, filter, unroll_inner, modulo_array, temp_array, modulo: true)
  if unroll_inner && !temp_array
    decl *tmp
  else
    decl tmp
  end
  index = lambda { |l = 0|
    indx = j + offset + k
    indx +=  l if unroll_inner
    if modulo
      if modulo_array
        indx = modulo_array[indx]
      else
        indx = Modulo(indx + n, n)
      end
    end
    indx
  }
  if unroll_inner
    unroll_inner.times { |l|
      pr tmp[l] === 0.0
    }
    pr For(k, 0, FILTER_SIZE-1, step: unroll_inner) {
      unroll_inner.times { |l|
        pr tmp[l] === tmp[l] + source[index.call(l), i] * filter[k + l]
      }
    }
    pr dest[i, j] === unroll_inner.times.collect { |l| tmp[l] }.reduce(:+)
  else
    pr tmp === 0.0
    pr For(k, 0, FILTER_SIZE-1) {
      pr tmp === tmp + source[index.call, i] * filter[k]
    }
    pr dest[i, j] === tmp
  end
end

def print_middle_loop(i, j, k, n, tmp, offset, source, dest, filter, unroll_inner, peel_middle, modulo_array, temp_array)
  if peel_middle
    pr For(j, 0, -offset - 1) {
      print_inner_loop(i, j, k, n, tmp, offset, source, dest, filter, unroll_inner, modulo_array, temp_array)
    }
    pr For(j, -offset, n - (FILTER_SIZE - 1 + offset) -1) {
      print_inner_loop(i, j, k, n, tmp, offset, source, dest, filter, unroll_inner, modulo_array, temp_array, modulo: false)
    }
    pr For(j, n - (FILTER_SIZE - 1 + offset), n - 1) {
      print_inner_loop(i, j, k, n, tmp, offset, source, dest, filter, unroll_inner, modulo_array, temp_array)
    }
  else
    pr For(j, 0, n-1) {
      print_inner_loop(i, j, k, n, tmp, offset, source, dest, filter, unroll_inner, modulo_array, temp_array)
    }
  end
end

def magic_filter_gen(t: false, unroll_inner: false, peel_middle: false, modulo_array: false, temp_array: true)
  raise "invalid unroll_inner: #{unroll_inner}" if unroll_inner && FILTER_SIZE % unroll_inner != 0
  n    = Int :n,    reference: true, dir: :in
  ndat = Int :ndat, reference: true, dir: :in
  source = Real :source, dim: [Dim(n), Dim(ndat)], dir: :in
  dest   = Real :dest,   dim: [Dim(ndat), Dim(n)], dir: :out
  if unroll_inner
    if temp_array
      tmp = Real :tmp, dim: Dim(unroll_inner), local: true
    else
      tmp = unroll_inner.times.collect { |i| Real :"tmp#{i}" }
    end
  else
    tmp = Real :tmp
  end
  i = Int :i
  j = Int :j
  k = Int :k
  offset = t ? -7 : -8
  if modulo_array
    modulo_array = Int :mod_arr, dim: Dim(offset, n + (FILTER_SIZE - 1 + offset) - 1), local: true
  end
  filter = Real :filt, dim: Dim(FILTER_SIZE), constant: ConstArray::new(t ? FILTER_REVERSE : FILTER), align: 16
  name = "magicfilter1d"
  name << "_t" if t
  name << "_ui#{unroll_inner}" if unroll_inner
  name << "_p" if peel_middle
  name << "_m" if modulo_array
  name << "_a" if temp_array
  magicfilter1d = Procedure(name, [n, ndat, source, dest]) {
    decl i, j, k
    if modulo_array
      decl modulo_array
      pr For(j, offset, n + (FILTER_SIZE - 1 + offset) - 1) {
        pr modulo_array[j] === Modulo(j + n, n)
      }
    end
    pr For(i, 0, ndat-1) {
      print_middle_loop(i, j, k, n, tmp, offset, source, dest, filter, unroll_inner, peel_middle, modulo_array, temp_array)
    }
  }
  k = CKernel::new {
    decl filter
    pr magicfilter1d
  }
  k.procedure = magicfilter1d
  k
end

def test_kernel(k, x, y, z, src, tmp, dest, ref)
  k.build(CFLAGS: CFLAGS)
  puts k.procedure.name
  apply_kernel(k, x, y, z, src, tmp, dest)
  err = (dest - ref).abs.max
  puts "Err: #{err}"
  raise "Wrong results!" if err > 1e-14
end


opt_space = OptimizationSpace::new( peel_middle: [false, true],
                                    modulo_array: [false, true],
                                    unroll_inner: [false, 1, 2, 4, 8, 16],
                                    temp_array: [true, false] )
optimizer = BruteForceOptimizer::new(opt_space, :randomize => true)


res = optimizer.optimize { |opt|
  puts opt
  k = magic_filter_gen(**opt)
  test_kernel(k, X, Y, Z, source, tmp, dest, dest_ref)
  r = bench_kernel(k, X, Y, Z, source, tmp, dest)
  puts "#{r} GFlop/s"
  - r
}
puts "Best: "
puts res[0]
puts magic_filter_gen(**res[0])
puts "#{-res[1]} GFlop/s"
