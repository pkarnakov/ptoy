#define CHECK_ERROR()                                                \
  do {                                                               \
    GLenum error = glGetError();                                     \
    if (error != GL_NO_ERROR) {                                      \
      std::string msg;                                               \
      do {                                                           \
        msg += reinterpret_cast<const char*>(gluErrorString(error)); \
        msg += "\n";                                                 \
        error = glGetError();                                        \
      } while (error != GL_NO_ERROR);                                \
      fassert(false, msg);                                           \
    }                                                                \
  } while (0);

std::string ReadFile(std::string path) {
  std::ifstream fs(path);
  fassert(fs.good(), "can't open file + '" + path + "'");
  std::stringstream ss;
  ss << fs.rdbuf();
  return ss.str();
}

namespace wrap {
auto glGetAttribLocation = [](GLuint prog, std::string varname,
                              bool nofail = false) {
  GLint loc = ::glGetAttribLocation(prog, varname.c_str());
  if (!nofail) {
    fassert(loc != -1, "attribute '" + varname + "' not found");
  }
  return loc;
};
auto glGetUniformLocation = [](GLuint prog, std::string varname) {
  GLint loc = ::glGetUniformLocation(prog, varname.c_str());
  fassert(loc != -1, "uniform '" + varname + "' not found");
  return loc;
};
void glCheckError() {
  GLenum e = glGetError();
  if (e != GL_NO_ERROR) {
    const std::string s(reinterpret_cast<const char*>(gluErrorString(e)));
    throw std::runtime_error(s);
  }
}

void glBufferDataReuse(
    GLenum target, GLsizeiptr size, const void* data, GLenum usage,
    GLuint vbo) {
  glBindBuffer(target, vbo);
  GLint cursize = 0;
  glGetBufferParameteriv(target, GL_BUFFER_SIZE, &cursize);

  if (cursize < size) { // FIXME: possible factor sizeof(GLfloat)
    glBufferData(target, size, data, usage);
  } else {
    glBufferSubData(target, 0, size, data);
  }
}

} // namespace wrap

template <class T>
GLenum GetGlTypeEnum();

template <>
GLenum GetGlTypeEnum<GLfloat>() {
  return GL_FLOAT;
}
template <>
GLenum GetGlTypeEnum<GLint>() {
  return GL_INT;
}
template <>
GLenum GetGlTypeEnum<GLchar>() {
  return GL_BYTE;
}

template <class T, int ncomp>
struct VertexAttribute {
  VertexAttribute(std::string name, GLuint program, bool nofail = true) {
    glGenBuffers(1, &id);
    location = wrap::glGetAttribLocation(program, name, nofail);
    if (location == -1) {
      return;
    }
    glEnableVertexAttribArray(location);
  }
  ~VertexAttribute() {
    glDeleteBuffers(1, &id);
  }
  void SetData(const std::vector<T>& data) const {
    if (location == -1) {
      return;
    }
    wrap::glBufferDataReuse(
        GL_ARRAY_BUFFER, data.size() * sizeof(T), data.data(), GL_DYNAMIC_DRAW,
        id);
    glVertexAttribPointer(
        location, ncomp, GetGlTypeEnum<T>(), GL_FALSE, ncomp * sizeof(T), NULL);
  }
  void SetData(const std::vector<std::array<T, ncomp>>& data) const {
    if (location == -1) {
      return;
    }
    wrap::glBufferDataReuse(
        GL_ARRAY_BUFFER, data.size() * sizeof(T) * ncomp, data.data(),
        GL_DYNAMIC_DRAW, id);
    glVertexAttribPointer(
        location, ncomp, GetGlTypeEnum<T>(), GL_FALSE, ncomp * sizeof(T), NULL);
  }
  void SetDataInt(const std::vector<T>& data) const {
    if (location == -1) {
      return;
    }
    wrap::glBufferDataReuse(
        GL_ARRAY_BUFFER, data.size() * sizeof(T), data.data(), GL_DYNAMIC_DRAW,
        id);
    glVertexAttribIPointer(
        location, ncomp, GetGlTypeEnum<T>(), ncomp * sizeof(T), NULL);
  }
  void SetDataInt(const std::vector<std::array<T, ncomp>>& data) const {
    if (location == -1) {
      return;
    }
    wrap::glBufferDataReuse(
        GL_ARRAY_BUFFER, data.size() * sizeof(T) * ncomp, data.data(),
        GL_DYNAMIC_DRAW, id);
    glVertexAttribIPointer(
        location, ncomp, GetGlTypeEnum<T>(), ncomp * sizeof(T), NULL);
  }
  GLuint id;
  GLint location = -1;
};

struct Texture {
  Texture(GLenum target_ = GL_TEXTURE_2D) : target(target_) {
    glGenTextures(1, &id);
    Bind();
    glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  }
  ~Texture() {
    glDeleteTextures(1, &id);
  }
  void SetData(
      const std::vector<float>& data, unsigned width, unsigned height) const {
    fassert_equal(data.size(), width * height);
    Bind();
    glTexImage2D(
        target, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, data.data());
  }
  void Bind() const {
    glBindTexture(target, id);
  }
  GLuint id;
  GLenum target;
};

struct Font {
  std::vector<float> mask;
  std::vector<float> xoffsets;
  std::vector<float> widths;
  size_t width;
  size_t height;

  Vect GetLength(const char* text) const {
    Scal res = 0;
    if (text) {
      while (*text) {
        res += widths[*text];
        ++text;
      }
    }
    return Vect(res, height);
  }
  Vect GetLength(const std::string& text) const {
    return GetLength(text.c_str());
  }
};

// Reads font mask and geometry from files generated by `font.py`.
Font LoadFont(std::string path) {
  Font font;

  {
    const std::string p = path + ".geom";
    std::ifstream f(p);
    fassert(f.good(), "Cannot open '" + p + "'");
    f >> font.width >> font.height;
    const size_t nchar = 256;
    font.xoffsets.resize(nchar);
    font.widths.resize(nchar);
    for (size_t ic = 0; ic < nchar; ++ic) {
      std::uint16_t xoffset, width;
      f >> xoffset >> width;
      font.xoffsets[ic] = xoffset;
      font.widths[ic] = width;
    }
  }

  {
    const std::string p = path + ".bin";
    std::ifstream f(p, std::ios::binary);
    fassert(f.good(), "Cannot open '" + p + "'");
    std::vector<unsigned char> buf(std::istreambuf_iterator<char>(f), {});
    fassert_equal(buf.size(), font.width * font.height);
    font.mask.resize(font.width * font.height);
    for (size_t i = 0; i < font.mask.size(); ++i) {
      font.mask[i] = 1 - buf[i] / 255.;
    }
  }
  return font;
}

void PrintShaderLog(GLuint shader) {
  if (glIsShader(shader)) {
    int max_length = 2048;
    int length = 0;

    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &max_length);

    char* log = new char[max_length];

    glGetShaderInfoLog(shader, max_length, &length, log);
    if (length > 0) {
      printf("%s\n", log);
    }

    delete[] log;
  } else {
    printf("Name %d is not a shader\n", shader);
  }
}

void PrintProgramLog(GLuint program) {
  if (glIsProgram(program)) {
    int max_length = 2048;
    int length = 0;

    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &max_length);

    char* log = new char[max_length];

    glGetProgramInfoLog(program, max_length, &length, log);
    if (length > 0) {
      printf("%s\n", log);
    }

    delete[] log;
  } else {
    printf("Name %d is not a program\n", program);
  }
}

GLuint CompileShaderSource(std::string src, GLenum type, std::string loc = "") {
  GLuint shader = glCreateShader(type);
  const GLchar* srcc[] = {src.c_str()};
  glShaderSource(shader, 1, srcc, NULL);
  glCompileShader(shader);

  GLint compiled = GL_FALSE;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
  if (compiled != GL_TRUE) {
    printf("Unable to compile shader %d in '%s'\n", shader, loc.c_str());
    PrintShaderLog(shader);
    fassert(false);
  }
  CHECK_ERROR();
  return shader;
}

GLuint CompileShaderFile(std::string path, GLenum type) {
  const std::string src = ReadFile(path);
  return CompileShaderSource(src, type, path);
}

GLuint CreateProgram(
    std::string vert_path, std::string frag_path, std::string geom_path) {
  GLuint program = glCreateProgram();

  std::vector<GLuint> todelete;

  if (vert_path.length()) {
    GLuint shader = CompileShaderFile(vert_path, GL_VERTEX_SHADER);
    glAttachShader(program, shader);
    todelete.push_back(shader);
    CHECK_ERROR();
  }

  if (frag_path.length()) {
    GLuint shader = CompileShaderFile(frag_path, GL_FRAGMENT_SHADER);
    glAttachShader(program, shader);
    todelete.push_back(shader);
    CHECK_ERROR();
  }

  if (geom_path.length()) {
    GLuint shader = CompileShaderFile(geom_path, GL_GEOMETRY_SHADER);
    glAttachShader(program, shader);
    todelete.push_back(shader);
    CHECK_ERROR();
  }

  glLinkProgram(program);
  { // error check
    GLint success = GL_TRUE;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (success != GL_TRUE) {
      PrintProgramLog(program);
      printf("Error linking program %d!\n", program);
      fassert(false);
    }
  }
  CHECK_ERROR();

  for (auto shader : todelete) {
    glDeleteShader(shader);
  }

  return program;
}
