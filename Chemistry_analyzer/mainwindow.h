#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProcess>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	MainWindow(QWidget* parent = nullptr);
	~MainWindow();

private slots:
	void on_pushButton_clicked();
	void on_action_H_triggered();
	void on_action_Qt_triggered();
	void on_action_code_triggered();
	void on_action_issue_triggered();
	void error();
	void output();
	void on_action_this_triggered();

private:
	Ui::MainWindow* ui;
	QProcess* core = new QProcess(this);
};
#endif // MAINWINDOW_H