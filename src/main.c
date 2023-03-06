#define CIMGUI_DEFINE_ENUMS_AND_STRUCTS
#define CIMGUI_USE_GLFW
#define CIMGUI_USE_OPENGL3
#include <cimgui.h>
#include <cimgui_impl.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <cglm/cglm.h>
// #define CIMGUI_DEFINE_ENUMS_AND_STRUCTS
// #include <cimgui/cimgui.h>

const char *vertexShaderSource =
		"#version 330 core\n"
		"layout (location = 0) in vec3 pos;\n"
		"void main()\n"
		"{\n"
		" gl_Position = vec4(pos.xyz, 1.0);\n"
		"}";
const char *initShaderSource =
		"#version 330 core\n"
		"uniform float iZoom;\n"
		"uniform vec2 iPosition;\n"
		"uniform vec2 iScreen;\n"
		"uniform vec2 iMouse;\n"
		"uniform vec2 iMove;\n"
		"uniform int iMode;\n"
		"in vec4 gl_FragCoord;\n"
		"out vec4 frag_color;\n"
		"#define E 2.71828182845904523536028747135266250\n"
		"#define ESCAPE 1000.\n"
		"#define PI 3.141592653\n"
		"vec2 cx_mul(vec2 a,vec2 b){\n"
		" return vec2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);\n"
		"}\n"
		"vec2 cx_sqr(vec2 a){\n"
		" float x2=a.x*a.x;\n"
		" float y2=a.y*a.y;\n"
		" float xy=a.x*a.y;\n"
		" return vec2(x2-y2,xy+xy);\n"
		"}\n"
		"vec2 cx_cube(vec2 a){\n"
		" float x2=a.x*a.x;\n"
		" float y2=a.y*a.y;\n"
		" float d=x2-y2;\n"
		" return vec2(a.x*(d-y2-y2),a.y*(x2+x2+d));\n"
		"}\n"
		"vec2 cx_div(vec2 a,vec2 b){\n"
		" float denom=1./(b.x*b.x+b.y*b.y);\n"
		" return vec2(a.x*b.x+a.y*b.y,a.y*b.x-a.x*b.y)*denom;\n"
		"}\n"
		"vec2 cx_sin(vec2 a){\n"
		" return vec2(sin(a.x)*cosh(a.y),cos(a.x)*sinh(a.y));\n"
		"}\n"
		"vec2 cx_cos(vec2 a){\n"
		" return vec2(cos(a.x)*cosh(a.y),-sin(a.x)*sinh(a.y));\n"
		"}\n"
		"vec2 cx_exp(vec2 a){\n"
		" return exp(a.x)*vec2(cos(a.y),sin(a.y));\n"
		"}\n"
		"vec2 fractal_f(vec2 z,vec2 c){\n"
		" return %s;\n"
		"}\n"
		"#if 1\n"
		"#define DO_LOOP(name)\\\n"
		"float smooth_i;\\\n"
		"for(i=0;i<ESCAPE;++i){\\\n"
		" vec2 ppz=pz;\\\n"
		" pz=z;\\\n"
		" z=name(z,c);\\\n"
		" if(dot(z,z)>ESCAPE){"
		"		const float mod = sqrt(dot(z, z));\\\n"
		"		smooth_i = float(i) - log2(max(1.0f, log2(mod)));\\\n"
		"		break;\\\n"
		" }\\\n"
		" ;\\\n"
		" "
		"}\n"
		"#else\n"
		"#define DO_LOOP(name)\\\n"
		"float smooth_i;\\\n"
		"for(i=0;i<ESCAPE;++i){\\\n"
		" z=name(z,c);\\\n"
		" if(dot(z,z)>ESCAPE){"
		"		const float mod = sqrt(dot(z, z));\\\n"
		"		smooth_i = float(i) - log2(max(1.0f, log2(mod)));\\\n"
		"		break;\\\n"
		" }\\\n"
		"}\n"
		"#endif\n"
		"vec3 gradient(float n){\n"
		" float div=1.f/ESCAPE;\n"
		" float red=10.f*n*div;\n"
		" float green=5.f*n*div-.5f;\n"
		" float blue=(6.f*n-9.f)/(2.f*(4.f*ESCAPE-6.f));\n"
		" return vec3(red,green,blue);\n"
		"}\n"
		"vec3 fractal(vec2 z,vec2 c){\n"
		" vec2 pz=z;\n"
		" int i;\n"
		" DO_LOOP(fractal_f);\n"
		" if (smooth_i) {\n"
		"  return gradient(smooth_i);"
		" }\n"
		" return gradient(i);\n"
		"}\n"
		"void main(){\n"
		" vec2 screen_pos=gl_FragCoord.xy-(iScreen.xy*.5);\n"
		" vec3 col=vec3(0.,0.,0.);\n"
		" vec2 c=vec2((screen_pos-iMove)*vec2(1.,-1.)*iZoom);\n"
		" if(iMode==0){\n"
		"   col+=fractal(c,c);\n"
		" }else if(iMode==1){\n"
		"   col+=fractal(c,c);\n"
		"   col+=fractal(vec2(screen_pos*vec2(1.,-1.)*.005434782608695652),iMouse);\n"
		" }else{\n"
		"   col+=fractal(c,iPosition);\n"
		" }\n"
		" if(iMode==1){\n"
		"   col*=.5;\n"
		" }\n"
		" gl_FragColor=vec4(clamp(col,0.,1.),1.);\n"
		"}";

int screen_width = 1000;
int screen_height = 700;

unsigned int VBO, VAO, EBO;
unsigned int shaderProgram;

double lastTime;
int nbFrames = 0;
int fps = 0;
char title[100];
char debugStr[100];

double whell = 0;
GLboolean rightClick = GL_FALSE;
GLboolean leftClick = GL_FALSE;

vec2 wherePressed = {0, 0};
bool canDrag = true;
vec2 mousePosition = {-100, -100};
vec2 cameraTranslation = {0, 0};
vec2 target = {0, 0};
double zoom = 184;

typedef enum
{
	mandelbrot = 0,
	view = 1,
	julia = 2
} modes;
modes mode = mandelbrot;

ImGuiIO *ioptr;

void gui_init(GLFWwindow *window)
{
	/* Initialize CIMGUI */
	// GL 3.2 + GLSL 130
	const char *glsl_version = "#version 130";

	igCreateContext(NULL);

	// set docking
	ioptr = igGetIO();
	ioptr->ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
// ioptr->ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;  // Enable Gamepad Controls
#ifdef IMGUI_HAS_DOCK
	ioptr->ConfigFlags |= ImGuiConfigFlags_DockingEnable;		// Enable Docking
	ioptr->ConfigFlags |= ImGuiConfigFlags_ViewportsEnable; // Enable Multi-Viewport / Platform Windows
#endif

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	igStyleColorsDark(NULL);
}

void gui_terminate(void)
{
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	igDestroyContext(NULL);
}

void gui_render(GLFWwindow *window)
{
	igRender();
	glfwMakeContextCurrent(window);
	ImGui_ImplOpenGL3_RenderDrawData(igGetDrawData());
}

#define BUF_SIZE 100
#define TOKEN_SIZE 4
#define AST_NEW(tag, ...) \
  ast_new((AST){tag, {.tag=(struct tag){__VA_ARGS__}}})
char input[BUF_SIZE] = "z^2+c";
char output[2 * BUF_SIZE] = "cx_sqr(z)+c";

typedef enum {f_sin, f_cos, f_exp} functions;
typedef enum {c_z, c_c, c_i, c_pi, c_tau, c_e} constants;

typedef struct AST AST; // Forward reference

struct AST
{
	enum
	{
		AST_IMAGINARY,
		AST_NUMBER,
		AST_FUNCTION,
		AST_CONSTANT,
		AST_UNR,
		AST_ADD,
		AST_SUB,
		AST_MUL,
		AST_DIV,
		AST_POW,
		AST_FUN,
	} tag;
	union
	{
		struct AST_IMAGINARY
		{
			float factor;
		} AST_IMAGINARY;
		struct AST_NUMBER
		{
			float number;
		} AST_NUMBER;
		struct AST_FUNCTION
		{
			functions function;
		} AST_FUNCTION;
		struct AST_CONSTANT
		{
			constants constant;
		} AST_CONSTANT;
		struct AST_UNR
		{
			AST *right;
		} AST_UNR;
		struct AST_ADD
		{
			AST *left;
			AST *right;
		} AST_ADD;
		struct AST_SUB
		{
			AST *left;
			AST *right;
		} AST_SUB;
		struct AST_MUL
		{
			AST *left;
			AST *right;
		} AST_MUL;
		struct AST_DIV
		{
			AST *left;
			AST *right;
		} AST_DIV;
		struct AST_POW
		{
			AST *left;
			AST *right;
		} AST_POW;
		struct AST_FUN
		{
			AST *left;
			AST *right;
		} AST_FUN;
	} data;
};

AST *ast_new(AST ast)
{
	AST *ptr = malloc(sizeof(AST));
	if (ptr)
		*ptr = ast;
	return ptr;
}

CGLM_INLINE
void print_float(float n, char *destStr)
{
	if (n == round(n))
	{
		snprintf(destStr, 2*BUF_SIZE, "%.1f", n);
	} else {
		snprintf(destStr, 2*BUF_SIZE, "%g", n);
	}
}

bool print_complex_mul(AST *left, AST *right)
{
	char floatStr[30];
	char *signStr = "%svec2(0.0, -%s)";
	char *unsignStr = "%svec2(0.0, %s)";
	char *str;
	AST *imag;
	AST *real;

	if (left->data.AST_CONSTANT.constant == c_i)
	{
		imag = left;
		real = right;
	} else if (right->data.AST_CONSTANT.constant == c_i) {
		imag = right;
		real = left;
	} else {
		return false;
	}

	if (real->tag == AST_UNR)
	{
		if (real->data.AST_UNR.right->tag != AST_NUMBER) return false;
		str = signStr;
		real = real->data.AST_UNR.right;
	}
	else
		str = unsignStr;

	if (real->tag != AST_NUMBER && real->tag != AST_UNR) return false;

	print_float(real->data.AST_NUMBER.number, floatStr);
	snprintf(output, 2*BUF_SIZE, str, output, floatStr);
	return true;
}

bool print_complex_add(AST *left, AST *right, bool sign)
{
	char floatStr1[30];
	char floatStr2[30];
	char *signStr1 = "%svec2(%s, -%s)";
	char *signStr2 = "%svec2(-%s, %s)";
	char *unsignStr = "%svec2(%s, %s)";
	char *str;
	AST *imag;
	AST *real;

	if (!sign)
		str = unsignStr;

	if (left->tag == AST_IMAGINARY || left->data.AST_CONSTANT.constant == c_i)
	{
		imag = left;
		real = right;
		if (sign) str = signStr2;
	} else if (right->tag == AST_IMAGINARY || right->data.AST_CONSTANT.constant == c_i) {
		imag = right;
		real = left;
		if (sign) str = signStr1;
	} else {
		return false;
	}

	if (real->tag != AST_NUMBER && real->tag != AST_UNR) return false;

	if (real->tag == AST_UNR)
	{
		if (real->data.AST_UNR.right->tag != AST_NUMBER) return false;
		str = signStr1;
		real = real->data.AST_UNR.right;
	}

	if (imag->tag == AST_CONSTANT)
	{
		print_float(real->data.AST_NUMBER.number, floatStr1);
		snprintf(output, 2*BUF_SIZE, str, output, floatStr1, "1.0");
		return true;
	}
	
	print_float(real->data.AST_NUMBER.number, floatStr1);
	print_float(imag->data.AST_IMAGINARY.factor, floatStr2);
	snprintf(output, 2*BUF_SIZE, str, output, floatStr1, floatStr2);
	return true;
}

void ast_print(AST *ptr) {
  AST ast = *ptr;
  switch (ast.tag) {
    case AST_IMAGINARY: {
      struct AST_IMAGINARY data = ast.data.AST_IMAGINARY;
			char floatStr[30];
			print_float(data.factor, floatStr);
			snprintf(output, 2*BUF_SIZE, "%svec2(0.0, %s)", output, floatStr);
      return;
    }
    case AST_NUMBER: {
      struct AST_NUMBER data = ast.data.AST_NUMBER;
			char floatStr[30];
			print_float(data.number, floatStr);
			snprintf(output, 2*BUF_SIZE, "%svec2(%s, 0.0)", output, floatStr);
      return;
    }
    case AST_FUNCTION: {
      struct AST_FUNCTION data = ast.data.AST_FUNCTION;

			switch (data.function)
			{
			case f_sin:
				snprintf(output, 2*BUF_SIZE, "%scx_sin", output);
				return;
			case f_cos:
				snprintf(output, 2*BUF_SIZE, "%scx_cos", output);
				return;
			case f_exp:
				snprintf(output, 2*BUF_SIZE, "%scx_exp", output);
				return;
			}
    }
    case AST_CONSTANT: {
      struct AST_CONSTANT data = ast.data.AST_CONSTANT;

			switch (data.constant)
			{
			case c_z:
				snprintf(output, 2*BUF_SIZE, "%sz", output);
				return;
			case c_c:
				snprintf(output, 2*BUF_SIZE, "%sc", output);
				return;
			case c_i:
				snprintf(output, 2*BUF_SIZE, "%svec2(0.0, 1.0)", output);
				return;
			case c_pi:
				snprintf(output, 2*BUF_SIZE, "%svec2(PI, 0.0)", output);
				return;
			case c_tau:
				snprintf(output, 2*BUF_SIZE, "%svec2(2.0 * PI, 0.0)", output);
				return;
			case c_e:
				snprintf(output, 2*BUF_SIZE, "%svec2(E, 0.0)", output);
				return;
			}
    }
    case AST_UNR: {
      struct AST_UNR data = ast.data.AST_UNR;
			snprintf(output, 2*BUF_SIZE, "%s-(", output);
      ast_print(data.right);
			snprintf(output, 2*BUF_SIZE, "%s)", output);
      return;
    }
    case AST_FUN: {
      struct AST_FUN data = ast.data.AST_FUN;
      ast_print(data.left);
			snprintf(output, 2*BUF_SIZE, "%s(", output);
      ast_print(data.right);
			snprintf(output, 2*BUF_SIZE, "%s)", output);
      return;
    }
    case AST_POW: {
      struct AST_POW data = ast.data.AST_POW;
			if (data.right->data.AST_NUMBER.number == 2)
			{
				snprintf(output, 2*BUF_SIZE, "%scx_sqr(", output);
				ast_print(data.left);
				snprintf(output, 2*BUF_SIZE, "%s)", output);
			} else if (data.right->data.AST_NUMBER.number == 3) {
				snprintf(output, 2*BUF_SIZE, "%scx_cube(", output);
				ast_print(data.left);
				snprintf(output, 2*BUF_SIZE, "%s)", output);
			}
      return;
    }
    case AST_ADD: {
      struct AST_ADD data = ast.data.AST_ADD;
			
			if (print_complex_add(data.left, data.right, false)) return;

			snprintf(output, 2*BUF_SIZE, "%s(", output);
      ast_print(data.left);
			snprintf(output, 2*BUF_SIZE, "%s + ", output);
      ast_print(data.right);
			snprintf(output, 2*BUF_SIZE, "%s)", output);
      return;
    }
    case AST_SUB: {
      struct AST_SUB data = ast.data.AST_SUB;
			
			if (print_complex_add(data.left, data.right, true)) return;

			snprintf(output, 2*BUF_SIZE, "%s(", output);
      ast_print(data.left);
			snprintf(output, 2*BUF_SIZE, "%s - ", output);
      ast_print(data.right);
			snprintf(output, 2*BUF_SIZE, "%s)", output);
      return;
    }
    case AST_MUL: {
      struct AST_MUL data = ast.data.AST_MUL;
			
			if (print_complex_mul(data.left, data.right)) return;

			snprintf(output, 2*BUF_SIZE, "%s(", output);
      ast_print(data.left);
			snprintf(output, 2*BUF_SIZE, "%s * ", output);
      ast_print(data.right);
			snprintf(output, 2*BUF_SIZE, "%s)", output);
      return;
    }
    case AST_DIV: {
      struct AST_DIV data = ast.data.AST_DIV;
			snprintf(output, 2*BUF_SIZE, "%scx_div(", output);
      ast_print(data.left);
			snprintf(output, 2*BUF_SIZE, "%s, ", output);
      ast_print(data.right);
			snprintf(output, 2*BUF_SIZE, "%s)", output);
      return;
    }
  }
}

// Top-level declarations
static AST *parse_add(void);
static const char *charPtr;

// Checks if state is a valid digit
static int is_letter(void) {
	return *charPtr >= 'a' && *charPtr <= 'z';
}

// Checks if state is a valid digit
static int is_digit(void) {
	return *charPtr >= '0' && *charPtr <= '9';
}

// Get the digit value of state
static int digit(void)
{
	return *charPtr - '0';
}

// Parses a number
static AST *number(void)
{
	float result = 0;

	while (is_digit())
	{
		int n = digit();

		result *= 10;
		result += n;

		++charPtr;
	}

	if (*charPtr == '.')
	{
		++charPtr;
		float dec = 0.1;

		while (is_digit())
		{
			int n = digit();

			result += n * dec;
			dec *= 0.1;

			++charPtr;
		}
	}

	if (*charPtr == 'i')
	{
		++charPtr; return AST_NEW(AST_IMAGINARY, result);
	} else
		return AST_NEW(AST_NUMBER, result);
}

// parses a constant
static AST *constant(char *conStr)
{
	constants con;

	if (!strcmp(conStr, "z"))
	{
		con = c_z;
	} else if (!strcmp(conStr, "c"))
	{
		con = c_c;
	} else if (!strcmp(conStr, "i"))
	{
		con = c_i;
	} else if (!strcmp(conStr, "pi"))
	{
		con = c_pi;
	} else if (!strcmp(conStr, "tau"))
	{
		con = c_tau;
	} else if (!strcmp(conStr, "e"))
	{
		con = c_e;
	} else {
		return parse_add();
	}

	return AST_NEW(AST_CONSTANT, con);
}

// parses a function
static AST *function(void)
{
	if (!is_letter()) return number();

	char funStr[TOKEN_SIZE];
	functions fun;
	int i = 0;

	while(is_letter())
	{
		if (i <= TOKEN_SIZE)
		{
			funStr[i] = *charPtr;
			funStr[i+1] = '\0';
		}

		++charPtr;
		++i;
	}

	if (!strcmp(funStr, "sin"))
	{
		fun = f_sin;
	} else if (!strcmp(funStr, "cos"))
	{
		fun = f_cos;
	} else if (!strcmp(funStr, "exp"))
	{
		fun = f_exp;
	} else {
		return constant(funStr);
	}

	AST *left = AST_NEW(AST_FUNCTION, fun);
	AST *right;

	if (*charPtr == '(')
	{
		++charPtr; // eat (
		right = parse_add();
		++charPtr; // eat )
		return AST_NEW(AST_FUN, left, right);
	}

	return parse_add();
}

// Parses a factor (unary, parenthesis, and number)
AST *parse_factor(void)
{
	switch (*charPtr)
	{
	case '+':
		++charPtr;
		return parse_factor();
	case '-':
		++charPtr;
		return AST_NEW(AST_UNR, parse_factor());
	case '(':
	{
		++charPtr; // eat (
		AST *result = parse_add();
		++charPtr; // eat )

		return result;
	}
	default:
		return function();
	}
}

// Parses ^
static AST *parse_expoent(void)
{
	AST *left = parse_factor();

	while (*charPtr == '^')
	{
		char op = *charPtr++;
		AST *right = parse_factor();

		left = AST_NEW(AST_POW, left, right);
	}

	return left;
}

// Parses * and /
static AST *parse_mul(void)
{
	AST *left = parse_expoent();

	while (*charPtr == '*' || *charPtr == '/')
	{
		char op = *charPtr++;
		AST *right = parse_expoent();

		if (op == '*')
			left = AST_NEW(AST_MUL, left, right);
		else
			left = AST_NEW(AST_DIV, left, right);
	}

	return left;
}

// Parses + and -
static AST *parse_add(void)
{
	AST *left = parse_mul();

	while (*charPtr == '+' || *charPtr == '-')
	{
		char op = *charPtr++;
		AST *right = parse_mul();

		if (op == '+')
			left = AST_NEW(AST_ADD, left, right);
		else
			left = AST_NEW(AST_SUB, left, right);
	}

	return left;
}

void reset(GLFWwindow *window)
{
	zoom = 184;
	glm_vec2_zero(cameraTranslation);
	glfwSetCursorPos(window, screen_width / 2, screen_height / 2);
}

char finalShaderSource[5000];

void compileShader(void)
{
	glDeleteProgram(shaderProgram);
	// vertex shader
	unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
	unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	snprintf(finalShaderSource, 5000, initShaderSource, output);
	const char* finalShaderSourceTmp = finalShaderSource;
	shaderProgram = glCreateProgram();
	glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
	glCompileShader(vertexShader);
	// fragment shader
	glShaderSource(fragmentShader, 1, &finalShaderSourceTmp, NULL);
	glCompileShader(fragmentShader);
	// link shaders
	glAttachShader(shaderProgram, vertexShader);
	glLinkProgram(shaderProgram);
	// link shaders
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);

	// mouse move
	glUseProgram(shaderProgram);
	int vertexLocation = glGetUniformLocation(shaderProgram, "iMove");
	glUniform2f(vertexLocation, 0, 0);
	vertexLocation = glGetUniformLocation(shaderProgram, "iZoom");
	glUniform1f(vertexLocation, 1 / zoom);
	vertexLocation = glGetUniformLocation(shaderProgram, "iMouse");
	glUniform2f(vertexLocation, 0, 0);
	vertexLocation = glGetUniformLocation(shaderProgram, "iScreen");
	glUniform2f(vertexLocation, screen_width, screen_height);
}

CGLM_INLINE
void convertStr(char *str, char *strDest)
{
	charPtr = str;
	AST *ast = parse_add();
	ast_print(ast);
}

int textInputCallback(ImGuiInputTextCallbackData* data) {
	memset(output, 0, 2 * BUF_SIZE);
	convertStr(data->Buf, output);

	return 0;
}

void gui_update(GLFWwindow *window)
{
	// start imgui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	igNewFrame();

	igBegin("Test", NULL, 0);
	igText("Fractal formula");
	igInputText("Formula", input, BUF_SIZE, ImGuiInputTextFlags_CallbackEdit, textInputCallback, NULL);
	igText("GLSL result: %s", output);
	if (igButton("compile", (ImVec2){100, 30}))
	{
		reset(window);
		compileShader();
	}

	igEnd();

	gui_render(window);

#ifdef IMGUI_HAS_DOCK
	if (ioptr->ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
	{
		GLFWwindow *backup_current_window = glfwGetCurrentContext();
		igUpdatePlatformWindows();
		igRenderPlatformWindowsDefault(NULL, NULL);
		glfwMakeContextCurrent(backup_current_window);
	}
#endif
}

CGLM_INLINE
void convertGrid(vec2 vector, vec2 dest)
{
	double x = (vector[0] + cameraTranslation[0] - screen_width / 2) / (zoom);
	double y = (vector[1] + cameraTranslation[1] - screen_height / 2) / (zoom);

	dest[0] = x;
	dest[1] = y;
}

void changeMode(modes toMod)
{
	// shaderProgram = shaders[toMod];
	// glUseProgram(shaderProgram);
	int vertexLocation = glGetUniformLocation(shaderProgram, "iMode");
	glUniform1i(vertexLocation, toMod);
	mode = toMod;
}

void resizeCallback(GLFWwindow *window, int width, int height)
{
	screen_width = width;
	screen_height = height;
	glViewport(0, 0, screen_width, screen_width);
	int vertexLocation = glGetUniformLocation(shaderProgram, "iScreen");
	glUniform2f(vertexLocation, screen_width, screen_height);
}

void scrollCallback(GLFWwindow *window, double xoffset, double yoffset)
{
	ImGuiIO *io = igGetIO();
	if (io->WantCaptureMouse)
		return;

	double toZoom = zoom * (0.5f * yoffset + 1);
	if (toZoom > 0)
	{
		vec2 vector1;
		vec2 vector2;
		glm_vec2_add(mousePosition, cameraTranslation, vector1);
		glm_vec2_sub(vector1, (vec2){screen_width / 2, screen_height / 2}, vector2);
		glm_vec2_copy(vector2, cameraTranslation);
		glfwSetCursorPos(window, screen_width / 2, screen_height / 2);
		glm_vec2_scale(cameraTranslation, toZoom / zoom, cameraTranslation);
		zoom = toZoom;
	}
}

void mouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
	ImGuiIO *io = igGetIO();
	if (io->WantCaptureMouse)
		return;

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		glm_vec2_copy(mousePosition, target);
		glm_vec2_copy(cameraTranslation, wherePressed);
		canDrag = true;
	}
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
	{
		if (mode == mandelbrot)
		{
			changeMode(view);
			int vertexLocation = glGetUniformLocation(shaderProgram, "iMove");
			glUniform2f(vertexLocation, cameraTranslation[0], cameraTranslation[1]);
			vertexLocation = glGetUniformLocation(shaderProgram, "iZoom");
			glUniform1f(vertexLocation, 1 / zoom);
		}
		else
		{
			changeMode(mandelbrot);
			reset(window);
		}
	}
	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE)
	{
		if (mode == view)
		{
			changeMode(julia);
			int vertexLocation = glGetUniformLocation(shaderProgram, "iPosition");
			vec2 position;
			convertGrid(mousePosition, position);
			glUniform2f(vertexLocation, position[0], position[1]);
			reset(window);
		}
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
	{
		canDrag = false;
	}
	if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS)
	{
		reset(window);
	}
}

void mouseCursorCallback(GLFWwindow *window, double xpos, double ypos)
{
	ImGuiIO *io = igGetIO();
	if (io->WantCaptureMouse)
		return;

	double x;
	double y;
	glfwGetCursorPos(window, &x, &y);
	glm_vec2_copy((vec2){x, y}, mousePosition);

	if (canDrag && glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) != GLFW_RELEASE)
	{
		vec2 draged;
		glm_vec2_sub(target, mousePosition, draged);
		glm_vec2_add(wherePressed, draged, cameraTranslation);
	}
}

void init()
{
	// fps counter
	lastTime = glfwGetTime();

	compileShader();

	// set up vertex data (and buffer(s)) and configure vertex attributes
	float vertices[] = {
			-1.0f, -1.0f, -0.0f,
			1.0f, 1.0f, -0.0f,
			-1.0f, 1.0f, -0.0f,
			1.0f, -1.0f, -0.0f};
	unsigned int indices[] = {
			0, 1, 2, // first Triangle
			0, 3, 1	 // second Triangle
	};

	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	// screen color
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
}

void render()
{
	double currentTime = glfwGetTime();

	glUseProgram(shaderProgram);
	glClear(GL_COLOR_BUFFER_BIT);
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

	// mouse move
	int vertexLocation = glGetUniformLocation(shaderProgram, "iMove");
	glUniform2f(vertexLocation, -cameraTranslation[0], cameraTranslation[1]);
	vertexLocation = glGetUniformLocation(shaderProgram, "iZoom");
	glUniform1f(vertexLocation, 1 / zoom);
	vertexLocation = glGetUniformLocation(shaderProgram, "iTime");
	glUniform1i(vertexLocation, floor(currentTime));
	if (mode == view)
	{
		vertexLocation = glGetUniformLocation(shaderProgram, "iMouse");
		vec2 position;
		convertGrid(mousePosition, position);
		glUniform2f(vertexLocation, position[0], position[1]);
	}

	// fps counter
	nbFrames++;
	if (currentTime - lastTime >= 1.0)
	{ // If last prinf() was more than 1 sec ago
		// printf and reset timer
		fps = nbFrames;
		nbFrames = 0;
		lastTime += 1.0;
	}

	// zoom += 0.01;
	// glm_vec2_scale(cameraTranslation, zoom/(zoom - 0.01), cameraTranslation);
}

int main(void)
{
	GLFWwindow *window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(screen_width, screen_height, "window", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);
	glfwSwapInterval(0);
	gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
	glfwSetCursorPosCallback(window, mouseCursorCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
	glfwSetScrollCallback(window, scrollCallback);
	glfwSetWindowSizeCallback(window, resizeCallback);

	printf("OpenGL loaded\nVersion: %s\n", glGetString(GL_VERSION));

	init();
	gui_init(window);
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		render();
		gui_update(window);
		sprintf(title, "fps: %d, %s", fps, output);
		glfwSetWindowTitle(window, title);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	glfwTerminate();
	gui_terminate();

	puts("fim do programa");

	return 0;
}