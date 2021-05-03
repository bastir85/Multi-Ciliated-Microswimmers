#ifndef ABBREV
#define ABBREV
// TODO move to cmake and gloabl header
#define TEST_CASES

// Defining global compiler abreviations...

#define PI 3.14159265358979323846
#ifdef DEBUG
#define PR(x) cout << #x "=" << x << endl;
#else
#define PR(x) ;
#endif

#define CP cout << "." << endl;
#define WRITE(x) out << #x << " " << x << endl;
#define REED(x)                                                                \
  in >> garb;                                                                  \
  if (garb != #x) {                                                            \
    cout << "error with " << garb << " as " << #x;                             \
    exit(0);                                                                   \
  }                                                                            \
  in >> x;
#define REEDS(x)                                                               \
  in >> garb;                                                                  \
  if (garb != #x) {                                                            \
    cout << "assume no " << #x << " defined. set it to 0! ";                   \
    x = 0;                                                                     \
  } else                                                                       \
    in >> x;
#define CHECK(x)                                                               \
  in >> garb;                                                                  \
  if (garb != #x) {                                                            \
    cout << "error with " << garb << " as " << #x;                             \
    exit(0);                                                                   \
  }                                                                            \
  in >> tempor;                                                                \
  if (tempor != x) {                                                           \
    cout << "error with" << #x;                                                \
    exit(0);                                                                   \
  }

#endif
