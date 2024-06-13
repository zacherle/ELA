#ifdef __cplusplus
extern"C" {
#endif

int filehyp_read(const char* hypname);
int filesta_read(const char* staname);
int filemod_read(const char* modname);

void setpar_model(const double model_err,const double reading_err, const char * modname);
void setpar_name_o(const char * hy3name);

#ifdef __cplusplus
}
#endif

