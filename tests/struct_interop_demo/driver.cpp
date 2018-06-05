#include <iostream>

extern "C" {
  struct Position { double x, y, z; };
  void set_position_array(Position *, const int);
}

int main() {
  const int n = 5;
  Position *p = new Position[n]; 

  set_position_array(p, n);

  for (int i=0; i<n; ++i)
    std::cout << p[i].x << p[i].y << p[i].z << std::endl;

  delete []p;
}
